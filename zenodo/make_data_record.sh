#!/usr/bin/env bash
# Builds the Zenodo data record folder from a completed run of run_analysis.py.
# Copies group-level result tables only and never touches the run outputs.
# Usage: bash zenodo/make_data_record.sh [run_output_dir] [record_dir]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${1:-${AMYG_PD_OUTPUT_DIR:?Set AMYG_PD_OUTPUT_DIR or pass the run output folder}}"
DST="${2:-$HOME/zenodo_amygdala_data_record}"
mkdir -p "$DST/tables" "$DST/figure_source_data"
python3 - "$SRC" "$DST" "$HERE" <<'PY'
import sys, os, re, shutil
import pandas as pd
src, dst, here = sys.argv[1:4]
# group-level tables written by analysis steps 02 to 11 (participant-level files are deliberately absent)
TABLES = ['participant_flow.csv','followup_duration_summary.csv','qc_outlier_group_breakdown.csv',
 'cs4_whole_amygdala_PD_vs_HC.csv','cs4_sectors_PD_vs_HC.csv','long5a_fixed_slopes_PD_vs_HC.csv',
 'long5b_whole_amygdala_interaction.csv','long5b_sectors_interaction.csv','long5b_nuclei_interaction_posthoc.csv',
 'demographics_table.csv','demographics_within_subtype_significance.csv','PD_vs_HC_clinical.csv',
 'VH11b_subtype_sectors.csv','VH11b_subtype_sectors_pairwise.csv','VH11b_subtype_nuclei_posthoc.csv',
 'VH11c_longitudinal_subtype_time.csv',
 'VH_9D_sectors.csv','VH_9D_nuclei.csv','VH_9E_within_group_slopes.csv','VH_9H_demographics.csv',
 'VH_9G_survival.csv','VH_9G_cohort_flow.csv',
 'CI_9D_sectors.csv','CI_9D_nuclei.csv','CI_9E_within_group_slopes.csv','CI_9H_demographics.csv',
 'CI_9G_survival.csv','CI_9G_cohort_flow.csv',
 'VHge2_9D_sectors.csv','VHge2_9D_nuclei.csv','VHge2_9G_survival.csv','VHge2_9G_cohort_flow.csv',
 'Step9J_random_slope_sensitivity.csv']
ID = re.compile(r'(^patno$|subject|participant_id|^id$|infodt|dob|date)', re.I)
DESC = {
 'region':'Amygdala region (whole amygdala, subregion or nucleus)','sector':'Amygdala subregion','col':'Internal column name of the region',
 'group':'Participant group','grp':'Participant group','subtype':'PD clinical subtype','analysis':'Analysis label',
 'n':'Number of participants','n_HC':'Number of healthy controls','n_PD':'Number of participants with PD',
 'n_early':'Number in Early group','n_neg':'Number in Non-Early group','n_pos':'Number in Early group','n_total':'Total participants',
 'mean_HC_mm3':'Mean volume in controls (mm3)','mean_PD_mm3':'Mean volume in PD (mm3)','HC_mean':'Mean in controls','PD_mean':'Mean in PD',
 't':'Welch t statistic','d':"Cohen's d (Early minus Non-Early, or PD minus HC)",'p':'Uncorrected p value','p_uncorrected':'Uncorrected p value',
 'p_fdr':'Benjamini-Hochberg FDR-corrected p value within the family defined in the article','KW_p':'Kruskal-Wallis p value',
 'KW_p_fdr':'Kruskal-Wallis p value, FDR-corrected','conover_p_fdr':'Conover post-hoc p value, FDR-corrected','comparison':'Pairwise comparison',
 'slope':'Annual volume change (mm3/year)','interaction':'Group-by-time interaction coefficient (mm3/year)','interaction_beta':'Group-by-time interaction coefficient (mm3/year)',
 'p_interaction':'p value of the group-by-time interaction','omnibus_chi2':'Likelihood-ratio chi-square for subtype-by-time (2 df)','p_omnibus':'p value of omnibus subtype-by-time test',
 'hr_per_SD_lower':'Hazard ratio per SD lower baseline volume','ci_lo':'Lower 95% confidence limit','ci_hi':'Upper 95% confidence limit',
 'logrank_p':'Log-rank p value with delayed entry','median_low':'Median years from diagnosis to onset, lower volume group','median_high':'Median years from diagnosis to onset, higher volume group',
 'p_ph':'Proportional hazards test p value','events':'Number of incident events','included':'Participants included in survival analysis',
 'with_valid_assessment':'Participants with a valid symptom assessment','prevalent_at_first':'Excluded, symptom present at first assessment','single_assessment':'Excluded, only one assessment',
 'tier':'Random-effects structure retained','converged':'Model converged','mean_years':'Mean follow-up (years)','sd_years':'SD of follow-up (years)',
 'median_years':'Median follow-up (years)','min_years':'Minimum follow-up (years)','max_years':'Maximum follow-up (years)',
 'n_with_followup':'Participants with at least two MRI visits','n_no_followup':'Participants with baseline MRI only','n_excluded':'Participants excluded','exclusion_type':'Outlier type',
 'Variable':'Clinical variable','variable':'Clinical variable','Early_mean_SD':'Early group mean and SD','NonEarly_mean_SD':'Non-Early group mean and SD',
 'symptom':'Symptom (VH or cognitive impairment)','year':'Years from first MRI','years_from_diagnosis':'Years from PD diagnosis',
 'emm_nonearly_mm3':'Estimated marginal mean, Non-Early group (mm3)','se_nonearly':'Standard error, Non-Early','emm_early_mm3':'Estimated marginal mean, Early group (mm3)','se_early':'Standard error, Early',
 'diff_early_minus_nonearly':'Difference in EMM, Early minus Non-Early (mm3)','survival':'Kaplan-Meier probability of remaining symptom-free','ci_lower':'Lower 95% confidence limit','ci_upper':'Upper 95% confidence limit',
 'n_at_risk':'Participants at risk','mean_mm3':'Mean volume (mm3)','sd_mm3':'SD of volume (mm3)','median_mm3':'Median volume (mm3)','q1_mm3':'First quartile (mm3)','q3_mm3':'Third quartile (mm3)',
}
rows, missing, blocked = [], [], []
def take(path, sub):
    df = pd.read_csv(path)
    bad = [c for c in df.columns if ID.search(str(c))]
    if bad: blocked.append((os.path.basename(path), bad)); return
    fn = os.path.basename(path); shutil.copy2(path, os.path.join(dst, sub, fn))
    for c in df.columns:
        rows.append({'file': f'{sub}/{fn}', 'column': c, 'type': str(df[c].dtype), 'description': DESC.get(str(c), '')})
for t in TABLES:
    p = os.path.join(src, t)
    take(p, 'tables') if os.path.exists(p) else missing.append(t)
fsd = os.path.join(src, 'figure_source_data')
if os.path.isdir(fsd):
    for f in sorted(os.listdir(fsd)):
        if f.endswith('.csv'): take(os.path.join(fsd, f), 'figure_source_data')
else: missing.append('figure_source_data/ (run step 13)')
xl = os.path.join(src, 'AMYGDALA_PD_MANUSCRIPT.xlsx')
if os.path.exists(xl): shutil.copy2(xl, dst)
pd.DataFrame(rows).to_csv(os.path.join(dst, 'data_dictionary.csv'), index=False)
for f in ['data_record_README.md']:
    shutil.copy2(os.path.join(here, f), os.path.join(dst, 'README.md'))
undesc = sum(1 for r in rows if not r['description'])
print(f'Copied {len({r["file"] for r in rows})} files to {dst}')
print(f'data_dictionary.csv: {len(rows)} columns, {undesc} still need a description')
if missing: print('Not found:', ', '.join(missing))
if blocked:
    print('BLOCKED (identifier-like columns, not copied):')
    for f, b in blocked: print(f'  {f}: {b}')
PY
echo "Done. Add descriptions for any blank entries in data_dictionary.csv before uploading to Zenodo."
