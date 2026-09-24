# Why Does Early PRBC Transfusion Improve Survival in Sepsis?

Analysis code for: *Why Does Early Packed Red Blood Cell Transfusion Improve Survival in Sepsis?
Mechanistic Mediation Analysis and Multicenter External Validation* (submitted to Journal of Advanced Research).

This repository contains the SQL extraction queries and Python analysis scripts used to produce all
results in the manuscript. **No patient-level data are included**: MIMIC-IV, eICU-CRD, and MIMIC-IV-Echo
are credentialed-access databases governed by PhysioNet data use agreements that prohibit redistribution.
Interested researchers can obtain access at https://physionet.org (CITI training + data use agreement required).

## Repository structure

```
sql/
  eicu_cohort_extraction.sql                     # eICU-CRD v2.0 cohort extraction (PostgreSQL)
  mimic_cohort_extraction_with_echo_lvef.sql     # MIMIC-IV cohort reconstruction + MIMIC-IV-Echo LVEF linkage
python/
  eicu_analysis.py                               # IPTW, mediation (with multiple imputation), subgroup analysis
  mimic_ef_subgroup_analysis.py                  # EF subgroup interaction analysis (ICD vs echocardiographic LVEF)
  competing_risk_analysis.py                     # Weighted Aalen-Johansen competing-risk analysis (7-day AKI, death as competing event)
```

## Reproduction pipeline

1. **Database setup**: Install PostgreSQL; load MIMIC-IV v3.1, eICU-CRD v2.0, and MIMIC-IV-Echo v1.0.1
   following the official build scripts from https://github.com/MIT-LCP/mimic-code
   and https://github.com/MIT-LCP/eicu-code (the MIMIC-IV derived concepts, `mimiciv_derived`,
   including Sepsis-3, first-day SOFA, and Charlson, must be built first).

2. **Cohort extraction**: Run the SQL scripts against the local databases to produce stay-level cohort
   tables with the variable definitions described in the manuscript Methods and Supplementary Table S1.

3. **Analysis**: Run the Python scripts on the extracted cohorts:
   ```
   python eicu_analysis.py <eicu_cohort.csv>
   python mimic_ef_subgroup_analysis.py <mimic_cohort.csv>
   ```

## Environment

- PostgreSQL 16+
- Python 3.12 with: pandas, numpy, statsmodels, scikit-learn, scipy

## Correspondence

【通讯作者邮箱】
