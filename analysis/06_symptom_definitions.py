"""Dated clinical visits, visual hallucination and cognitive impairment status, early versus non-early groups."""

_updrs3 = safe_read_csv(f'{BASE_DIR}/PPMI_Clinical_13Apr2026/Motor___MDS-UPDRS/MDS-UPDRS_Part_III_12Apr2026.csv', low_memory=False)
_updrs3.columns = [c.replace('-','_') for c in _updrs3.columns]
_updrs3['PATNO'] = _updrs3['PATNO'].astype(str)
_updrs3['INFODT'] = pd.to_datetime(_updrs3['INFODT'], errors='coerce')
_date_lookup = _updrs3[['PATNO','EVENT_ID','INFODT']].dropna(subset=['INFODT']).drop_duplicates(['PATNO','EVENT_ID'])
_date_lookup['Visit'] = _date_lookup['EVENT_ID'].map(EVENT_MAP)

_img_bl = df.sort_values('years_from_first_mri').drop_duplicates('PATNO', keep='first')[['PATNO','Visit','years_from_first_mri']]
_img_bl_dated = _img_bl.merge(_date_lookup[['PATNO','Visit','INFODT']], on=['PATNO','Visit'], how='left')
_n_matched = _img_bl_dated['INFODT'].notna().sum()
print(f"Imaging baseline visits matched to a clinical INFODT date: {_n_matched} / {len(_img_bl_dated)}")
_first_mri_date_lookup = _img_bl_dated.dropna(subset=['INFODT']).set_index('PATNO')['INFODT'].to_dict()

_date_lookup['first_mri_date'] = _date_lookup['PATNO'].map(_first_mri_date_lookup)
_date_lookup['years_from_first_mri_CLINICAL'] = (
    (_date_lookup['INFODT'] - _date_lookup['first_mri_date']).dt.days / 365.25
)
print(f"Clinical visits with a computed years_from_first_mri_CLINICAL: "
      f"{_date_lookup['years_from_first_mri_CLINICAL'].notna().sum()} / {len(_date_lookup)}")

_clin_cols_dated = ['U1_NP1HALL','HVLT_DVT_TOTAL_RECALL','HVLT_DVT_RECOG_DISC_INDEX',
                    'SEMFLU_DVS_SFTANIM','SDMT_DVT_SDM','LNS_DVS_LNS']
_clin_cols_dated_present = [c for c in _clin_cols_dated if c in mega.columns]
_mega_sub_dated = mega[['PATNO','EVENT_ID'] + _clin_cols_dated_present].drop_duplicates(['PATNO','EVENT_ID'])
_final_dated = _date_lookup.merge(_mega_sub_dated, on=['PATNO','EVENT_ID'], how='left')

_final_dated_clean = _final_dated[['PATNO','EVENT_ID','Visit','INFODT','first_mri_date',
                                   'years_from_first_mri_CLINICAL'] + _clin_cols_dated_present
                                  ].sort_values(['PATNO','years_from_first_mri_CLINICAL'])
_final_dated_clean.to_excel(f'{OUT_DIR}/clinical_visit_times_all_visits.xlsx', index=False)
print(f"Saved: {OUT_DIR}/clinical_visit_times_all_visits.xlsx "
      f"({len(_final_dated_clean)} rows, {_final_dated_clean['PATNO'].nunique()} subjects)")

_followup_dated = _final_dated_clean[_final_dated_clean['years_from_first_mri_CLINICAL'] >= 0].copy()
_followup_dated.to_excel(f'{OUT_DIR}/clinical_visit_times_followup_only.xlsx', index=False)
print(f"Saved: {OUT_DIR}/clinical_visit_times_followup_only.xlsx "
      f"({len(_followup_dated)} rows, {_followup_dated['PATNO'].nunique()} subjects)")

_dx = safe_read_csv(f'{BASE_DIR}/PPMI_Clinical_13Apr2026/Medical_History/PD_Diagnosis_History_23Jun2026.csv', low_memory=False)
_dx['PATNO'] = _dx['PATNO'].astype(str)
_dx['PDDXDT'] = pd.to_datetime(_dx['PDDXDT'], errors='coerce')
_dx_lookup = _dx[['PATNO','PDDXDT']].dropna(subset=['PDDXDT']).drop_duplicates('PATNO')
print(f"Subjects with a valid PDDXDT: {len(_dx_lookup)}")

_first_mri_lookup = pd.read_excel(f'{OUT_DIR}/clinical_visit_times_all_visits.xlsx')[['PATNO','first_mri_date']].dropna().drop_duplicates('PATNO')
_first_mri_lookup['PATNO'] = _first_mri_lookup['PATNO'].astype(str)
_first_mri_lookup['first_mri_date'] = pd.to_datetime(_first_mri_lookup['first_mri_date'])

_dd_new = _first_mri_lookup.merge(_dx_lookup, on='PATNO', how='inner')
_dd_new['disease_duration_at_first_mri'] = (_dd_new['first_mri_date'] - _dd_new['PDDXDT']).dt.days / 365.25
print(f"Subjects with both first_mri_date and PDDXDT: {len(_dd_new)}")

_all_visits_dx = pd.read_excel(f'{OUT_DIR}/clinical_visit_times_all_visits.xlsx')
_all_visits_dx['PATNO'] = _all_visits_dx['PATNO'].astype(str)
_merged_dx = _all_visits_dx.merge(_dd_new[['PATNO','disease_duration_at_first_mri']], on='PATNO', how='inner')
_merged_dx['years_from_diagnosis'] = _merged_dx['years_from_first_mri_CLINICAL'] + _merged_dx['disease_duration_at_first_mri']
_merged_dx_fu = _merged_dx[_merged_dx['years_from_first_mri_CLINICAL'] >= 0].copy()
print(f"Visits (follow-up only): {len(_merged_dx_fu)}, subjects: {_merged_dx_fu['PATNO'].nunique()}")

def _build_early_groups(data, outcome_col, label):
    avail = data.dropna(subset=[outcome_col]).copy()
    pos = avail[avail[outcome_col]==1].copy()
    first_onset = pos.sort_values('years_from_diagnosis').drop_duplicates('PATNO', keep='first')
    first_onset = first_onset[['PATNO','years_from_diagnosis']].rename(columns={'years_from_diagnosis':'onset_years_from_dx'})
    last_fu = avail.groupby('PATNO')['years_from_diagnosis'].max().reset_index()
    last_fu.columns = ['PATNO','last_dx_followup_years']
    subj = avail[['PATNO']].drop_duplicates().merge(first_onset, on='PATNO', how='left').merge(last_fu, on='PATNO', how='left')
    subj['group'] = 'excluded_insufficient_followup'
    subj.loc[(subj['onset_years_from_dx'].notna()) & (subj['onset_years_from_dx'] <= 5), 'group'] = 'early_converter'
    subj.loc[(subj['onset_years_from_dx'].notna()) & (subj['onset_years_from_dx'] > 5), 'group'] = 'not_early_late_converter'
    subj.loc[(subj['onset_years_from_dx'].isna()) & (subj['last_dx_followup_years'] >= 5), 'group'] = 'not_early_confirmed_free'
    subj['early'] = np.where(subj['group']=='early_converter', 1,
                     np.where(subj['group'].isin(['not_early_late_converter','not_early_confirmed_free']), 0, np.nan))
    print(f"\n{label} groups:")
    print(subj['group'].value_counts())
    print(f"Final: early={int((subj['early']==1).sum())}, not_early={int((subj['early']==0).sum())}, excluded={subj['early'].isna().sum()}")
    return subj

_vh_dx = _merged_dx_fu.dropna(subset=['U1_NP1HALL']).copy()
_vh_dx['U1_NP1HALL'] = pd.to_numeric(_vh_dx['U1_NP1HALL'], errors='coerce')
print(f"NP1HALL = 101 (not applicable) visits dropped from VH groups: {int((_vh_dx['U1_NP1HALL']==101).sum())}")
_vh_dx = _vh_dx[_vh_dx['U1_NP1HALL'].notna() & (_vh_dx['U1_NP1HALL'] != 101)].copy()
_vh_dx['VH_pos'] = (_vh_dx['U1_NP1HALL'] >= 1).astype(float)
vh_early_groups = _build_early_groups(_vh_dx.rename(columns={'VH_pos':'outcome_pos'}), 'outcome_pos', 'VH (NP1HALL>=1)')
vh_early_groups.to_excel(f'{OUT_DIR}/VH_early_groups_FINAL.xlsx', index=False)

if 'DVS_JLO_MSSA' not in _merged_dx_fu.columns:
    _j9c = safe_read_csv(f'{BASE_DIR}/PPMI_Clinical_13Apr2026/Non-motor_Assessments/Benton_Judgement_of_Line_Orientation_12Apr2026.csv', low_memory=False)
    _j9c.columns = [c.replace('-', '_') for c in _j9c.columns]
    _j9c['DVS_JLO_MSSA'] = pd.to_numeric(_j9c['DVS_JLO_MSSA'], errors='coerce')
    _pdt = _merged_dx_fu['PATNO'].dtype
    _j9c['PATNO'] = _j9c['PATNO'].astype(str) if _pdt == object else _j9c['PATNO'].astype(_pdt)
    _j9c = _j9c.dropna(subset=['DVS_JLO_MSSA']).drop_duplicates(['PATNO', 'EVENT_ID'])[['PATNO', 'EVENT_ID', 'DVS_JLO_MSSA']]
    _merged_dx_fu = _merged_dx_fu.merge(_j9c, on=['PATNO', 'EVENT_ID'], how='left')
    print(f"  JLO merged onto _merged_dx_fu (9C): {int(_merged_dx_fu['DVS_JLO_MSSA'].notna().sum())} visits with a score")
_ci_clin_cols = ['HVLT_DVT_TOTAL_RECALL','HVLT_DVT_RECOG_DISC_INDEX','SEMFLU_DVS_SFTANIM','SDMT_DVT_SDM','LNS_DVS_LNS','DVS_JLO_MSSA']
_ci_clin_present = [c for c in _ci_clin_cols if c in _merged_dx_fu.columns]
print("  9C CI tests used:", _ci_clin_present)
_ci_dx = _merged_dx_fu.dropna(subset=_ci_clin_present, how='all').copy()
for _nm9, _c9, _th9 in [('imp_HVLT_recall','HVLT_DVT_TOTAL_RECALL',35), ('imp_HVLT_recog','HVLT_DVT_RECOG_DISC_INDEX',35),
                        ('imp_SemFluency','SEMFLU_DVS_SFTANIM',5.5), ('imp_SDMT','SDMT_DVT_SDM',35),
                        ('imp_LNS','LNS_DVS_LNS',5.5), ('imp_JLO','DVS_JLO_MSSA',5.5)]:
    if _c9 in _ci_dx.columns:
        _x9 = pd.to_numeric(_ci_dx[_c9], errors='coerce')
        _ci_dx[_nm9] = np.where(_x9.notna(), (_x9 <= _th9).astype(float), np.nan)
_imp_cols_dx = [c for c in ['imp_HVLT_recall','imp_HVLT_recog','imp_SemFluency','imp_SDMT','imp_LNS','imp_JLO'] if c in _ci_dx.columns]
_ci_dx['n_avail'] = _ci_dx[_imp_cols_dx].notna().sum(axis=1)
_ci_dx['n_impaired'] = _ci_dx[_imp_cols_dx].sum(axis=1)
_ci_dx['CI_pos'] = np.where(_ci_dx['n_avail'] >= 4, (_ci_dx['n_impaired'] >= 2).astype(float), np.nan)
ci_early_groups = _build_early_groups(_ci_dx.dropna(subset=['CI_pos']).rename(columns={'CI_pos':'outcome_pos'}), 'outcome_pos', 'CI')
ci_early_groups.to_excel(f'{OUT_DIR}/CI_early_groups_FINAL.xlsx', index=False)
