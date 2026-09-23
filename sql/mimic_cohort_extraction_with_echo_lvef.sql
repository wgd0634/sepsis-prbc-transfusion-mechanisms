-- =====================================================================
-- MIMIC-IV 脓毒症输血队列重建 + 超声 LVEF 链接（本地 mimiciv 库）
-- 队列定义严格按论文 Methods：
--   成人 + 首次 ICU + ICU>24h + Sepsis-3(mimiciv_derived.sepsis3)
--   基线 = suspected_infection_time
--   排除: 基线前死亡 / 基线前 PRBC / 基线前 AKI(Cr>=1.5x 前7天最低值)
--   保留: 基线后24h内最低 Hb <= 10
--   暴露: 基线后24h内 PRBC (inputevents itemid 225168/226368/220996)
--   结局: 28天/90天死亡 (patients.dod), 7天AKI (KDIGO Cr, 参照=基线后24h最低Cr)
--   EF: ICD (4282/I502=HFrEF, 4283/I503=HFpEF) + 超声 LVEF
-- =====================================================================
SET search_path TO mimiciv_derived, mimiciv_hosp, mimiciv_icu;

CREATE TEMP TABLE echo_studies_tmp (subject_id BIGINT, study_id BIGINT, study_datetime TIMESTAMP);
\COPY echo_studies_tmp FROM 'C:/Users/wgd06/Documents/kimi/tasks/2026-09-19/22-43-19-87e50a1c/pancreatitis/echo_studies.csv' DELIMITER ',' CSV HEADER;

CREATE TEMP TABLE lvef_tmp (measurement_id TEXT, lvef TEXT);
\COPY lvef_tmp FROM 'C:/Users/wgd06/Documents/kimi/tasks/2026-09-19/22-43-19-87e50a1c/pancreatitis/lvef_values.csv' DELIMITER ',' CSV HEADER;

-- 每项研究的 LVEF（measurement_id 关联）
CREATE TEMP TABLE study_lvef AS
SELECT e.subject_id, e.study_id, e.study_datetime,
       CAST(v.lvef AS numeric) AS lvef
FROM echo_studies_tmp e
JOIN lvef_tmp v ON v.measurement_id = CAST(e.measurement_id AS TEXT)
WHERE v.lvef ~ '^\d+\.?\d*$';

WITH base AS (
    SELECT s.stay_id, s.subject_id, d.hadm_id, d.gender, d.admission_age AS age,
           d.dod, d.icu_intime, d.icu_outtime, d.los_icu,
           s.suspected_infection_time AS baseline_t
    FROM mimiciv_derived.sepsis3 s
    JOIN mimiciv_derived.icustay_detail d ON d.stay_id = s.stay_id
    WHERE d.first_icu_stay IS TRUE
      AND d.admission_age >= 18
      AND d.los_icu > 1
      AND NOT (d.dod IS NOT NULL AND d.dod < s.suspected_infection_time)   -- 基线前死亡
),
prbc_all AS (
    SELECT b.stay_id,
           MIN(CASE WHEN ie.starttime < b.baseline_t THEN ie.starttime END) AS prbc_before_t,
           MIN(CASE WHEN ie.starttime >= b.baseline_t
                     AND ie.starttime < b.baseline_t + interval '24 hours'
                    THEN ie.starttime END) AS prbc_in24h_t,
           SUM(CASE WHEN ie.starttime >= b.baseline_t
                     AND ie.starttime < b.baseline_t + interval '24 hours'
                    THEN ie.amount ELSE 0 END) AS prbc_ml_24h
    FROM base b
    JOIN inputevents ie ON ie.stay_id = b.stay_id
     AND ie.itemid IN (225168, 226368, 220996)     -- PRBC
    GROUP BY b.stay_id
),
aki_pre AS (
    -- 基线前7天 Cr：KDIGO 1.5x 最低值规则
    SELECT b.stay_id,
           MAX(le.valuenum) AS scr_max_pre,
           MIN(le.valuenum) AS scr_min_pre
    FROM base b
    JOIN labevents le ON le.subject_id = b.subject_id
     AND le.itemid = 50912
     AND le.charttime BETWEEN b.baseline_t - interval '7 days' AND b.baseline_t
     AND le.valuenum > 0
    GROUP BY b.stay_id
),
hb AS (
    SELECT b.stay_id, MIN(le.valuenum) AS hb_min_24h
    FROM base b
    JOIN labevents le ON le.subject_id = b.subject_id
     AND le.itemid = 51222                       -- Hemoglobin
     AND le.valuenum > 0
     AND le.charttime BETWEEN b.baseline_t AND b.baseline_t + interval '24 hours'
    GROUP BY b.stay_id
),
cohort AS (
    SELECT b.*, hb.hb_min_24h,
           CASE WHEN p.prbc_in24h_t IS NOT NULL THEN 1 ELSE 0 END AS prbc_24h_flag,
           p.prbc_ml_24h
    FROM base b
    JOIN hb ON hb.stay_id = b.stay_id
    LEFT JOIN prbc_all p ON p.stay_id = b.stay_id
    LEFT JOIN aki_pre ak ON ak.stay_id = b.stay_id
    WHERE hb.hb_min_24h <= 10
      AND p.prbc_before_t IS NULL
      AND NOT COALESCE(ak.scr_max_pre >= 1.5 * ak.scr_min_pre, false)
),
-- ===== 协变量 =====
sofa AS (
    SELECT c.stay_id, f.sofa AS sofa_score
    FROM cohort c
    LEFT JOIN first_day_sofa f ON f.stay_id = c.stay_id
),
cx AS (
    SELECT c.stay_id,
           MAX(ch.charlson_comorbidity_index) AS charlson,
           MAX(CASE WHEN ch.renal_disease = 1 THEN 1 ELSE 0 END) AS ckd,
           MAX(CASE WHEN ch.diabetes_without_cc = 1 OR ch.diabetes_with_cc = 1 THEN 1 ELSE 0 END) AS diabetes
    FROM cohort c
    LEFT JOIN charlson ch ON ch.hadm_id = c.hadm_id
    GROUP BY c.stay_id
),
vent AS (
    SELECT DISTINCT c.stay_id, 1 AS mech_vent_24h
    FROM cohort c
    JOIN ventilation v ON v.stay_id = c.stay_id
     AND v.ventilation_status = 'InvasiveVent'
     AND v.starttime < c.baseline_t + interval '24 hours'
     AND v.endtime   > c.baseline_t
),
vaso AS (
    SELECT DISTINCT c.stay_id, 1 AS vasopressor_24h
    FROM cohort c
    JOIN LATERAL (
        SELECT 1 FROM norepinephrine   x WHERE x.stay_id = c.stay_id AND x.starttime < c.baseline_t + interval '24 hours' AND x.endtime > c.baseline_t
        UNION ALL SELECT 1 FROM epinephrine     x WHERE x.stay_id = c.stay_id AND x.starttime < c.baseline_t + interval '24 hours' AND x.endtime > c.baseline_t
        UNION ALL SELECT 1 FROM vasopressin     x WHERE x.stay_id = c.stay_id AND x.starttime < c.baseline_t + interval '24 hours' AND x.endtime > c.baseline_t
        UNION ALL SELECT 1 FROM dopamine        x WHERE x.stay_id = c.stay_id AND x.starttime < c.baseline_t + interval '24 hours' AND x.endtime > c.baseline_t
        UNION ALL SELECT 1 FROM phenylephrine   x WHERE x.stay_id = c.stay_id AND x.starttime < c.baseline_t + interval '24 hours' AND x.endtime > c.baseline_t
    ) x ON true
),
rrt24 AS (
    SELECT DISTINCT c.stay_id, 1 AS rrt_24h
    FROM cohort c
    LEFT JOIN rrt r ON r.stay_id = c.stay_id
     AND r.charttime BETWEEN c.baseline_t AND c.baseline_t + interval '24 hours'
     AND r.dialysis_active = 1
    LEFT JOIN crrt k ON k.stay_id = c.stay_id
     AND k.charttime BETWEEN c.baseline_t AND c.baseline_t + interval '24 hours'
    WHERE r.stay_id IS NOT NULL OR k.stay_id IS NOT NULL
),
-- ===== EF 定义 1: ICD 编码 =====
icd_hf AS (
    SELECT c.stay_id,
           MAX(CASE WHEN (dg.icd_version = 9 AND dg.icd_code LIKE '4282%')
                     OR (dg.icd_version = 10 AND dg.icd_code LIKE 'I502%') THEN 1 ELSE 0 END) AS sys_hf,
           MAX(CASE WHEN (dg.icd_version = 9 AND dg.icd_code LIKE '4283%')
                     OR (dg.icd_version = 10 AND dg.icd_code LIKE 'I503%') THEN 1 ELSE 0 END) AS dia_hf,
           MAX(CASE WHEN (dg.icd_version = 9 AND dg.icd_code LIKE '428%')
                     OR (dg.icd_version = 10 AND dg.icd_code LIKE 'I50%') THEN 1 ELSE 0 END) AS any_hf
    FROM cohort c
    JOIN diagnoses_icd dg ON dg.hadm_id = c.hadm_id
    GROUP BY c.stay_id
),
-- ===== EF 定义 2: 超声 LVEF（基线-7d ~ +24h 内离基线最近的一次）====
echo_near AS (
    SELECT DISTINCT ON (c.stay_id)
        c.stay_id, v.lvef,
        v.study_datetime,
        ABS(EXTRACT(EPOCH FROM (v.study_datetime - c.baseline_t))) AS dist_sec
    FROM cohort c
    JOIN study_lvef v ON v.subject_id = c.subject_id
     AND v.study_datetime BETWEEN c.baseline_t - interval '7 days' AND c.baseline_t + interval '24 hours'
    ORDER BY c.stay_id, dist_sec
),
-- ===== 结局 =====
aki_out AS (
    SELECT c.stay_id,
           MIN(CASE WHEN le.charttime <= c.baseline_t + interval '24 hours' THEN le.valuenum END) AS scr_baseline,
           MAX(le.valuenum) AS scr_max_7d
    FROM cohort c
    JOIN labevents le ON le.subject_id = c.subject_id
     AND le.itemid = 50912
     AND le.valuenum > 0
     AND le.charttime BETWEEN c.baseline_t AND c.baseline_t + interval '7 days'
    GROUP BY c.stay_id
)
SELECT c.stay_id, c.subject_id, c.gender, c.age,
       c.hb_min_24h, sf.sofa_score AS sofa, cx.charlson,
       COALESCE(cx.ckd,0) AS ckd, COALESCE(cx.diabetes,0) AS diabetes,
       COALESCE(vt.mech_vent_24h,0) AS mech_vent_24h,
       COALESCE(vs.vasopressor_24h,0) AS vasopressor_24h,
       COALESCE(rt.rrt_24h,0) AS rrt_24h,
       c.prbc_24h_flag, c.prbc_ml_24h,
       e.lvef, e.study_datetime AS echo_datetime,
       COALESCE(ic.sys_hf,0) AS sys_hf, COALESCE(ic.dia_hf,0) AS dia_hf, COALESCE(ic.any_hf,0) AS any_hf,
       CASE WHEN c.dod IS NOT NULL AND c.dod <= c.baseline_t + interval '28 days' THEN 1 ELSE 0 END AS death_28d,
       CASE WHEN c.dod IS NOT NULL AND c.dod <= c.baseline_t + interval '90 days' THEN 1 ELSE 0 END AS death_90d,
       CASE WHEN ao.scr_max_7d - ao.scr_baseline >= 0.3
              OR ao.scr_max_7d >= 1.5 * ao.scr_baseline THEN 1 ELSE 0 END AS aki_7d
FROM cohort c
LEFT JOIN sofa sf ON sf.stay_id = c.stay_id
LEFT JOIN cx   cx ON cx.stay_id = c.stay_id
LEFT JOIN vent vt ON vt.stay_id = c.stay_id
LEFT JOIN vaso vs ON vs.stay_id = c.stay_id
LEFT JOIN rrt24 rt ON rt.stay_id = c.stay_id
LEFT JOIN icd_hf ic ON ic.stay_id = c.stay_id
LEFT JOIN echo_near e ON e.stay_id = c.stay_id
LEFT JOIN aki_out ao ON ao.stay_id = c.stay_id
ORDER BY c.stay_id;
