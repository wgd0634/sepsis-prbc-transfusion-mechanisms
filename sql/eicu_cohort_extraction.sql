-- =====================================================================
-- eICU-CRD 2.0 (本地 PostgreSQL, public schema) 脓毒症输血验证队列提取
-- 基于 eICU提取代码.sql 适配：
--   1) eicu_crd. -> public（无前缀）
--   2) icd9code 为逗号分隔多值编码，正则改为按编码 token 匹配（去除点/空格后按逗号边界）
--   3) PRBC celllabel 补 '%prbc%'/'%red blood cell%'（本地主标签为 pRBCs/PRBC）
--   4) cellvaluenumber -> cellvaluenumeric；cellpath 实际以 'flowsheet|...' 开头，
--      入量/出量判定改为 LIKE '%intake (ml)%' / '%output (ml)%'
--   5) drugrate 为 varchar，CAST 前用正则守卫；NE 剂量纳入 levophed，排除 Volume 行
--   6) apacheversion 本地为 IV/IVa 两种，均纳入
--   7) age 数值转换加正则守卫
-- =====================================================================
WITH first_icu AS (
    SELECT p.patientunitstayid AS stay_id, p.uniquepid, p.gender, p.age,
           p.hospitaldischargestatus, p.hospitaldischargeoffset,
           p.unitdischargeoffset
    FROM patient p
    WHERE p.unitvisitnumber = 1                       -- 首次ICU
      AND p.unitdischargeoffset > 1440                -- 停留>24h
),
sepsis AS (
    -- Angus ICD 方案 + 诊断字符串兜底（icd9code 多值逗号分隔，按 token 匹配）
    SELECT d.patientunitstayid AS stay_id,
           MIN(d.diagnosisoffset) AS sepsis_offset    -- 分钟
    FROM diagnosis d
    WHERE regexp_replace(lower(d.icd9code), '[\.\s]', '', 'g')
          ~ '(^|,)(038\d*|99591|99592|78552|r652\d*|a020|a021|a027|a220|a221|a227|a267|a327|a40\d*|a41\d*|b377\d*)'
       OR lower(d.diagnosisstring) LIKE '%sepsis%'
       OR lower(d.diagnosisstring) LIKE '%septicemia%'
       OR lower(d.diagnosisstring) LIKE '%septic shock%'
    GROUP BY d.patientunitstayid
),
cohort AS (
    SELECT f.*, s.sepsis_offset,
           CASE WHEN f.age = '> 89' THEN 90
                WHEN f.age ~ '^\d+(\.\d+)?$' THEN CAST(f.age AS numeric)
           END AS age_num
    FROM first_icu f
    JOIN sepsis s ON f.stay_id = s.stay_id
    WHERE s.sepsis_offset IS NOT NULL AND s.sepsis_offset >= 0
      AND f.age <> '' AND (f.age = '> 89' OR (f.age ~ '^\d+(\.\d+)?$' AND CAST(f.age AS numeric) >= 18))
),
hb AS (
    SELECT c.stay_id, MIN(l.labresult) AS hb_min_24h
    FROM cohort c
    JOIN lab l ON l.patientunitstayid = c.stay_id
     AND l.labname = 'Hgb'
     AND l.labresult > 0
     AND l.labresultoffset BETWEEN c.sepsis_offset AND c.sepsis_offset + 1440
    GROUP BY c.stay_id
),
prbc AS (
    SELECT c.stay_id,
           MIN(CASE WHEN io.intakeoutputoffset < c.sepsis_offset
                    THEN io.intakeoutputoffset END) AS prbc_before_offset,
           MIN(CASE WHEN io.intakeoutputoffset >= c.sepsis_offset
                     AND io.intakeoutputoffset < c.sepsis_offset + 1440
                    THEN io.intakeoutputoffset END) AS prbc_in24h_offset,
           SUM(CASE WHEN io.intakeoutputoffset >= c.sepsis_offset
                     AND io.intakeoutputoffset < c.sepsis_offset + 1440
                    THEN io.cellvaluenumeric ELSE 0 END) AS prbc_ml_24h
    FROM cohort c
    JOIN intakeoutput io ON io.patientunitstayid = c.stay_id
     AND (lower(io.celllabel) LIKE '%prbc%'
          OR lower(io.celllabel) LIKE '%packed red blood%'
          OR lower(io.celllabel) LIKE '%red blood cell%')
    GROUP BY c.stay_id
),
aki_before AS (
    -- 诊断前窗口 = 入ICU至诊断时点（eICU无院前数据）
    SELECT c.stay_id,
           MAX(l.labresult) AS scr_max_pre,
           MIN(l.labresult) AS scr_min_pre
    FROM cohort c
    JOIN lab l ON l.patientunitstayid = c.stay_id
     AND l.labname = 'creatinine'
     AND l.labresultoffset BETWEEN 0 AND c.sepsis_offset
    GROUP BY c.stay_id
),
base_cohort AS (
    SELECT c.stay_id, c.gender, c.age_num AS age, c.sepsis_offset,
           c.hospitaldischargestatus,
           hb.hb_min_24h,
           CASE WHEN p.prbc_in24h_offset IS NOT NULL THEN 1 ELSE 0 END AS prbc_24h_flag,
           p.prbc_ml_24h
    FROM cohort c
    JOIN hb              hb ON c.stay_id = hb.stay_id
    LEFT JOIN prbc        p ON c.stay_id = p.stay_id
    LEFT JOIN aki_before ak ON c.stay_id = ak.stay_id
    WHERE hb.hb_min_24h <= 10
      AND p.prbc_before_offset IS NULL
      AND NOT COALESCE(ak.scr_max_pre >= 1.5 * ak.scr_min_pre, false)
),
-- ===== 协变量 =====
apache AS (
    SELECT patientunitstayid AS stay_id, apachescore
    FROM apachepatientresult
    WHERE apacheversion IN ('IV', 'IVa')
),
comorb AS (
    SELECT b.stay_id,
           MAX(CASE WHEN regexp_replace(lower(d.icd9code), '[\.\s]', '', 'g') ~ '(^|,)(585|n18)' THEN 1 ELSE 0 END) AS ckd,
           MAX(CASE WHEN regexp_replace(lower(d.icd9code), '[\.\s]', '', 'g') ~ '(^|,)(250|e1[01])' THEN 1 ELSE 0 END) AS diabetes,
           MAX(CASE WHEN regexp_replace(lower(d.icd9code), '[\.\s]', '', 'g') ~ '(^|,)(4282|i502)' THEN 1 ELSE 0 END) AS sys_hf,
           MAX(CASE WHEN regexp_replace(lower(d.icd9code), '[\.\s]', '', 'g') ~ '(^|,)(4283|i503)' THEN 1 ELSE 0 END) AS dia_hf,
           MAX(CASE WHEN regexp_replace(lower(d.icd9code), '[\.\s]', '', 'g') ~ '(^|,)(428|i50)'  THEN 1 ELSE 0 END) AS any_hf
    FROM base_cohort b
    JOIN diagnosis d ON d.patientunitstayid = b.stay_id
    GROUP BY b.stay_id
),
vent AS (
    SELECT b.stay_id, 1 AS mech_vent_24h
    FROM base_cohort b
    JOIN treatment t ON t.patientunitstayid = b.stay_id
     AND lower(t.treatmentstring) LIKE '%ventilat%'
     AND t.treatmentoffset < b.sepsis_offset + 1440
    GROUP BY b.stay_id
),
vaso AS (
    SELECT b.stay_id, 1 AS vasopressor_24h
    FROM base_cohort b
    JOIN infusiondrug id ON id.patientunitstayid = b.stay_id
     AND (lower(id.drugname) LIKE '%norepinephrine%' OR lower(id.drugname) LIKE '%levophed%'
          OR lower(id.drugname) LIKE '%epinephrine%'
          OR lower(id.drugname) LIKE '%vasopressin%' OR lower(id.drugname) LIKE '%dopamine%'
          OR lower(id.drugname) LIKE '%phenylephrine%')
     AND id.infusionoffset < b.sepsis_offset + 1440
    GROUP BY b.stay_id
),
rrt AS (
    SELECT b.stay_id, 1 AS rrt_24h
    FROM base_cohort b
    JOIN treatment t ON t.patientunitstayid = b.stay_id
     AND (lower(t.treatmentstring) LIKE '%dialysis%' OR lower(t.treatmentstring) LIKE '%crrt%'
          OR lower(t.treatmentstring) LIKE '%hemofiltr%')
     AND t.treatmentoffset < b.sepsis_offset + 1440
    GROUP BY b.stay_id
),
-- ===== 中介变量 =====
fluids AS (
    SELECT b.stay_id,
           SUM(CASE WHEN lower(io.cellpath) LIKE '%intake (ml)%' THEN io.cellvaluenumeric ELSE 0 END) AS total_input_ml_24h,
           SUM(CASE WHEN lower(io.cellpath) LIKE '%output (ml)%' THEN io.cellvaluenumeric ELSE 0 END) AS total_output_ml_24h
    FROM base_cohort b
    JOIN intakeoutput io ON io.patientunitstayid = b.stay_id
     AND io.intakeoutputoffset >= b.sepsis_offset
     AND io.intakeoutputoffset <  b.sepsis_offset + 1440
    GROUP BY b.stay_id
),
lac0 AS (
    SELECT DISTINCT ON (b.stay_id) b.stay_id, l.labresult AS lactate_0h
    FROM base_cohort b
    JOIN lab l ON l.patientunitstayid = b.stay_id
     AND l.labname = 'lactate'
     AND l.labresultoffset BETWEEN b.sepsis_offset - 360 AND b.sepsis_offset + 360
    ORDER BY b.stay_id, ABS(l.labresultoffset - b.sepsis_offset)
),
lac24 AS (
    SELECT DISTINCT ON (b.stay_id) b.stay_id, l.labresult AS lactate_24h
    FROM base_cohort b
    JOIN lab l ON l.patientunitstayid = b.stay_id
     AND l.labname = 'lactate'
     AND l.labresultoffset BETWEEN b.sepsis_offset + 720 AND b.sepsis_offset + 2160
    ORDER BY b.stay_id, ABS(l.labresultoffset - (b.sepsis_offset + 1440))
),
ne0 AS (
    SELECT b.stay_id,
           MAX(CASE WHEN id.drugrate ~ '^\d+\.?\d*$' THEN CAST(id.drugrate AS numeric) END) AS norepi_0h
    FROM base_cohort b
    JOIN infusiondrug id ON id.patientunitstayid = b.stay_id
     AND (lower(id.drugname) LIKE '%norepinephrine%' OR lower(id.drugname) LIKE '%levophed%')
     AND lower(id.drugname) NOT LIKE '%volume%'
     AND id.infusionoffset >= b.sepsis_offset
     AND id.infusionoffset <  b.sepsis_offset + 360
    GROUP BY b.stay_id
),
ne24 AS (
    SELECT b.stay_id,
           MAX(CASE WHEN id.drugrate ~ '^\d+\.?\d*$' THEN CAST(id.drugrate AS numeric) END) AS norepi_24h
    FROM base_cohort b
    JOIN infusiondrug id ON id.patientunitstayid = b.stay_id
     AND (lower(id.drugname) LIKE '%norepinephrine%' OR lower(id.drugname) LIKE '%levophed%')
     AND lower(id.drugname) NOT LIKE '%volume%'
     AND id.infusionoffset >= b.sepsis_offset + 720
     AND id.infusionoffset <  b.sepsis_offset + 2160
    GROUP BY b.stay_id
),
-- ===== 结局 =====
aki_out AS (
    SELECT b.stay_id,
           MIN(CASE WHEN l.labresultoffset <= b.sepsis_offset + 1440
                    THEN l.labresult END) AS scr_baseline,
           MAX(l.labresult) AS scr_max_7d
    FROM base_cohort b
    JOIN lab l ON l.patientunitstayid = b.stay_id
     AND l.labname = 'creatinine'
     AND l.labresultoffset BETWEEN b.sepsis_offset AND b.sepsis_offset + 10080
    GROUP BY b.stay_id
)
SELECT b.stay_id, b.gender, b.age, b.sepsis_offset AS sepsis_offset_min,
       b.hb_min_24h, ap.apachescore,
       COALESCE(cm.ckd,0) AS ckd, COALESCE(cm.diabetes,0) AS diabetes,
       COALESCE(cm.sys_hf,0) AS sys_hf, COALESCE(cm.dia_hf,0) AS dia_hf, COALESCE(cm.any_hf,0) AS any_hf,
       COALESCE(vt.mech_vent_24h,0) AS mech_vent_24h,
       COALESCE(vs.vasopressor_24h,0) AS vasopressor_24h,
       COALESCE(rt.rrt_24h,0) AS rrt_24h,
       b.prbc_24h_flag, b.prbc_ml_24h,
       f.total_input_ml_24h, f.total_output_ml_24h,
       COALESCE(f.total_input_ml_24h,0) - COALESCE(f.total_output_ml_24h,0) AS fluid_balance_ml_24h,
       l0.lactate_0h, l24.lactate_24h, n0.norepi_0h, n24.norepi_24h,
       CASE WHEN b.hospitaldischargestatus = 'Expired' THEN 1 ELSE 0 END AS death_inhosp,
       CASE WHEN ao.scr_max_7d - ao.scr_baseline >= 0.3
              OR ao.scr_max_7d >= 1.5 * ao.scr_baseline THEN 1 ELSE 0 END AS aki_7d
FROM base_cohort b
LEFT JOIN apache ap ON b.stay_id = ap.stay_id
LEFT JOIN comorb cm ON b.stay_id = cm.stay_id
LEFT JOIN vent   vt ON b.stay_id = vt.stay_id
LEFT JOIN vaso   vs ON b.stay_id = vs.stay_id
LEFT JOIN rrt    rt ON b.stay_id = rt.stay_id
LEFT JOIN fluids  f ON b.stay_id = f.stay_id
LEFT JOIN lac0   l0 ON b.stay_id = l0.stay_id
LEFT JOIN lac24 l24 ON b.stay_id = l24.stay_id
LEFT JOIN ne0    n0 ON b.stay_id = n0.stay_id
LEFT JOIN ne24  n24 ON b.stay_id = n24.stay_id
LEFT JOIN aki_out ao ON b.stay_id = ao.stay_id
ORDER BY b.stay_id;
