# Source code

Analysis scripts: recipient age and post-transplant outcomes, donor-derived cell-free DNA (dd-cfDNA), and de novo
donor-specific antibodies (dnDSA).

Scripts expect input files in `data/` and write results to `output/`.

## Structure

| Folder | Script | Description |
|---|---|---|
| `R/` | `helpers.R` | Shared plotting theme, palettes and export functions |
| `data_preparation/` | `extract_donor_info.py` | Extracts donor variables from organ offers |
| | `calculate_hla_mismatch.R` | HLA mismatch level (1-4) from the A:B:DR mismatch string |
| | `calculate_risk_indices.R` | UK Kidney Donor and Recipient Risk Index and quartiles |
| `age_outcomes/` | `baseline_and_outcomes.R` | Baseline characteristics, clinical outcomes (OR, RR), survival; age <60 vs ≥60 |
| | `egfr_mixed_models.R` | eGFR over the first year, linear mixed-effects models |
| | `length_of_stay.R` | Length of stay, negative binomial regression |
| `ddcfdna/` | `cohort_comparison.R` | Baseline comparison of the two dd-cfDNA cohorts |
| | `result_categories.R` | dd-cfDNA results by clinical category |
| | `sample_quality.R` | Purity and concentration by collection tube and timepoint |
| `dnDSA/` | `firth_regression.R` | Univariable and multivariable Firth logistic regression, dnDSA trajectories |
| | `ML_prespecified_models.R` | Pre-specified predictors from a previous study; Firth regression and XGBoost |

## Input files

| File | Used by |
|---|---|
| `data/donor_reports/*.pdf` | `extract_donor_info.py` |
| `data/eligibility.xlsx` | `calculate_hla_mismatch.R` |
| `data/donors.xlsx` | `calculate_risk_indices.R` |
| `data/recipients.csv`, `data/recipient_donor_indices.csv` | `calculate_risk_indices.R` |
| `data/cohort.xlsx` | `baseline_and_outcomes.R`, `length_of_stay.R` |
| `data/cohort_long.xlsx` | `egfr_mixed_models.R` (one row per patient and timepoint, first sheet) |
| `data/ddcfdna.csv` | `ddcfdna/` scripts (first row is a title row) |
| `data/ddcfdna_qc_exclusions.csv` | `result_categories.R`, `sample_quality.R` (columns `cf_id`, `timepoint`, `analysis`) |
| `data/dsa_cohort.xlsx` | `dnDSA/` scripts |

## Running

Run all scripts from the repository root, for example:

```
Rscript data_preparation/calculate_hla_mismatch.R
Rscript age_outcomes/egfr_mixed_models.R          # full cohort
Rscript age_outcomes/egfr_mixed_models.R 70-79    # over 70s recipient age subgroup analysis
python data_preparation/extract_donor_info.py data/donor_reports output/donors_extracted.xlsx
```

Missing R packages are installed from CRAN on first run. The main packages are
lme4, lmerTest, logistf, MASS, xgboost, pROC, survival, survminer, tableone,
gtsummary, flextable, officer, rvg, ggplot2 and transplantr. The donor information
extraction script needs Python 3 with pandas, openpyxl and pdfplumber.

Tables are written as Word documents and figures as PowerPoint slides with
editable vector graphics, with PNG or PDF copies where noted in each script.
