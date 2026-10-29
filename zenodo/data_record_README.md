# Source data: Amygdala subregion atrophy in Parkinson's disease is associated with visual hallucinations and cognitive impairment

**Authors** Ajay Konuri, Jack Anderson, Lachlan Churchill, Anna Ignatavicius, Glenda M Halliday, Simon JG Lewis, Ramon Landin-Romero, Elie Matar
**Contact** Ajay Konuri, Brain and Mind Centre, University of Sydney ([email])
**Version** 1.0.0 · **Licence** CC BY 4.0
**Code** https://doi.org/10.5281/zenodo.[CODE_RECORD] · **Preprint** [medRxiv DOI] · **Article** [journal DOI]

## Description
Group-level results behind every table and figure in the article. The study examined whole amygdala, subregion (sCLR, CMA, BLA) and nucleus volumes in 468 participants with Parkinson's disease and 165 healthy controls from the Parkinson's Progression Markers Initiative (PPMI), and their associations with early and incident visual hallucinations, cognitive impairment and clinical subtypes.

This record contains summary statistics and model outputs only. It contains no participant-level values, identifiers or dates.

## Contents
| Path | Content |
|---|---|
| `AMYGDALA_PD_MANUSCRIPT.xlsx` | All tables below as one workbook, with an INDEX sheet |
| `tables/` | One CSV per result table (participant flow, PD vs HC, Early vs Non-Early VH and cognitive impairment, survival, subtypes, sensitivity analyses) |
| `figure_source_data/` | Values plotted in Figures 1 to 3 (baseline group summaries, estimated marginal means by year, Kaplan-Meier estimates with numbers at risk) |
| `data_dictionary.csv` | Every column in every file, with type and description |

## Table to manuscript map
| File | Manuscript item |
|---|---|
| tables/participant_flow.csv, followup_duration_summary.csv, qc_outlier_group_breakdown.csv | Methods, Table 1 |
| tables/demographics_table.csv, demographics_within_subtype_significance.csv, PD_vs_HC_clinical.csv | Table 1 |
| tables/cs4_*.csv, long5a_*.csv, long5b_*.csv | Results, PD versus HC |
| tables/VH_9*.csv | Results, Early VH; Supplementary Table 1; Figures 1 and 3A |
| tables/CI_9*.csv | Results, early cognitive impairment; Supplementary Table 2; Figures 2 and 3B |
| tables/VHge2_9*.csv | Results, stricter VH threshold |
| tables/Step9J_random_slope_sensitivity.csv | Results, random-slope sensitivity |
| tables/VH11b_*.csv, VH11c_*.csv | Results, clinical subtypes; Table 2 |
| figure_source_data/fig1*, fig2*, fig3* | Figures 1, 2 and 3 |

## Conventions
Volumes are in mm³ after adjustment for total intracranial volume, age and sex. Slopes and interaction coefficients are in mm³ per year. Hazard ratios are per SD lower baseline volume. `p_fdr` columns are Benjamini-Hochberg corrected within the families defined in the article. Empty cells are missing values.

## Access to participant-level data
Participant-level clinical and imaging data cannot be redistributed under the PPMI Data Use Agreement. Researchers can apply for access at https://www.ppmi-info.org/access-data-specimens/download-data (RRID:SCR_006431). Imaging data were downloaded in June 2026 and clinical data in April 2026. Running the archived code on a PPMI download reproduces every file here.

## Funding
This research was funded in whole or in part by Aligning Science Across Parkinson's (ASAP) (https://ror.org/03zj4c476) [grant number ASAP-027449] through the Michael J. Fox Foundation for Parkinson's Research (MJFF) (https://ror.org/03arq3225). For the purpose of open access, the authors have applied a CC BY 4.0 public copyright license to all Author Accepted Manuscripts arising from this submission.
