# Amygdala subregion atrophy in Parkinson's disease

Analysis code for

Konuri A, Anderson J, Churchill L, Ignatavicius A, Halliday GM, Lewis SJG, Landin-Romero R, Matar E. *Amygdala subregion atrophy in Parkinson's disease is associated with visual hallucinations and cognitive impairment.* npj Parkinson's Disease (under review).

The code takes harmonised FreeSurfer amygdala nucleus volumes and clinical data from the Parkinson's Progression Markers Initiative (PPMI) and reproduces every statistical result, table and figure in the manuscript.

## Contents

| File | Manuscript item |
|---|---|
| `analysis/common.py` | Paths, constants and shared helpers |
| `analysis/01_load_qc_adjust.py` | Outlier detection and adjustment for intracranial volume, age and sex (Methods) |
| `analysis/02_cohort_summary.py` | Participant flow and follow-up duration (Table 1) |
| `analysis/03_amygdala_subregions.py` | Grouping of nuclei into sCLR, CMA and BLA (Methods) |
| `analysis/04_pd_vs_hc.py` | PD versus controls, cross-sectional and longitudinal (Results) |
| `analysis/05_clinical_tables.py` | Demographic and clinical characteristics (Table 1) |
| `analysis/06_symptom_definitions.py` | Visual hallucination and cognitive impairment status, early and non-early groups (Methods) |
| `analysis/07_early_symptom_groups.py` | Early versus non-early VH and CI (Results, Supplementary Tables 1 and 2) |
| `analysis/08_incident_survival.py` | Cox and Kaplan-Meier analyses of incident VH and CI, sensitivity analysis restricting VH to formed hallucinations, MDS-UPDRS item 1.2 score of 2 or more (Results) |
| `analysis/09_figures.py` | Figures 1 to 3 |
| `analysis/10_random_slope_sensitivity.py` | Random-slope sensitivity models |
| `analysis/11_clinical_subtypes.py` | Clinical subtype comparisons (Results, Table 2) |
| `analysis/12_export_workbook.py` | All result tables in one Excel workbook |
| `analysis/13_export_figure_source_data.py` | Group-level source data for Figures 1 to 3 (Zenodo data record) |
| `run_analysis.py` | Runs the steps in order |
| `KeyResourcesTable.csv` | ASAP Key Resources Table |
| `zenodo/` | Kit for building the Zenodo data record |

The steps share variables and must be run through `run_analysis.py`.

## Requirements

- Python 3.12 or later (tested with 3.14.7)
- Packages in `requirements.txt` (pandas, NumPy, SciPy, statsmodels, lifelines, scikit-learn, scikit-posthocs, matplotlib, openpyxl)
- Runs on macOS and Linux. A full run takes a few minutes on a standard laptop.

## Installation

```bash
git clone https://github.com/akonuri/amygdala-subregions-pd.git
cd amygdala-subregions-pd
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Data

No participant data are included in this repository.

Data were obtained from the PPMI database (https://www.ppmi-info.org/access-data-specimens/download-data, RRID:SCR_006431), downloaded in June 2026. Researchers can access PPMI data after applying through the PPMI website and agreeing to the PPMI Data Use Agreement.

The analysis expects the following files in one folder.

| File | Contents |
|---|---|
| `harmonised_with_corrected_subtypes.csv` | FreeSurfer amygdala nucleus, cortical and subcortical volumes per participant and visit, harmonised across scanner manufacturers with longCombat, with diagnostic group, age, sex, visit, years from first MRI and clinical subtype |
| `PPMI_volumes_clinical.csv` | Estimated total intracranial volume and clinical scale totals per participant and visit |
| `PPMI_mega_clinical.csv` | Merged PPMI clinical and neuropsychological data per participant and visit |
| `PPMI_Clinical_13Apr2026/` | PPMI study data download (MDS-UPDRS Part III, Benton Judgement of Line Orientation, socio-economics, PD diagnosis history) |

Amygdala nuclei were segmented with FreeSurfer `recon-all` and `segmentHA_T1.sh`.

Group-level source data for all figures and tables are available on Zenodo (https://doi.org/10.5281/zenodo.[DATA_RECORD]). Build that record from a full run with `bash zenodo/make_data_record.sh` (see `zenodo/ZENODO_GUIDE.md`).

**Upstream steps not in this repository.** [Add or link the FreeSurfer `recon-all`/`segmentHA_T1.sh` batch script, the longCombat R script and the clinical subtype classification script that produce `harmonised_with_corrected_subtypes.csv`.]

## Usage

```bash
export AMYG_PD_DATA_DIR=/path/to/input/folder
export AMYG_PD_OUTPUT_DIR=/path/to/output/folder   # optional, defaults to a subfolder of the input folder
python run_analysis.py 2>&1 | tee analysis_log.txt
```

The log reports every statistic in the manuscript. Result tables are written as CSV and Excel files, figures to `publication_plots/`, figure source data to `figure_source_data/`, and all tables are combined in `AMYGDALA_PD_MANUSCRIPT.xlsx`.

## Methods summary

Cross-sectional outliers were identified at the first MRI using Mahalanobis distance in a principal component space fitted to control data, and longitudinal outliers from subject-specific slopes more than four robust median absolute deviations from the group median. Retained volumes were adjusted for intracranial volume, age and sex. Group comparisons used Welch's t-tests and Kruskal-Wallis tests with Conover post-hoc tests. Longitudinal change was modelled with linear mixed-effects models with participant random intercepts. Incident symptoms were modelled with Cox proportional hazards models with delayed entry at the first clinical assessment. False discovery rate correction was applied within each family of tests. Full details are in the manuscript.

## Citation

If you use this code, please cite the article above and the archived software release (Zenodo DOI [CODE DOI]).

## Licence

MIT. See `LICENSE`.

## Acknowledgements

This research was funded in whole or in part by Aligning Science Across Parkinson's (ASAP) (https://ror.org/03zj4c476) [grant number ASAP-027449] through the Michael J. Fox Foundation for Parkinson's Research (MJFF) (https://ror.org/03arq3225). For the purpose of open access, the authors have applied a CC BY 4.0 public copyright license to all Author Accepted Manuscripts arising from this submission. Data used in the preparation of this work were obtained from the Parkinson's Progression Markers Initiative (PPMI) database. PPMI is a public-private partnership funded by the Michael J. Fox Foundation for Parkinson's Research and funding partners.

## Contact

Ajay Konuri, Brain and Mind Centre, University of Sydney. Questions about the code can be raised as GitHub issues.
