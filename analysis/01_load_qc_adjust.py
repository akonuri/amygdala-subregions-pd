"""Load harmonised volumes and clinical data, detect cross-sectional and longitudinal outliers, adjust volumes for eTIV, age and sex."""

harm = safe_read_csv(f'{BASE_DIR}/harmonised_with_corrected_subtypes.csv')
harm.columns = [c.replace('-','_') for c in harm.columns]
harm['PATNO'] = harm['PATNO'].astype(str)

vol = safe_read_csv(f'{BASE_DIR}/PPMI_volumes_clinical.csv', low_memory=False)
vol.columns = [c.replace('-','_') for c in vol.columns]
vol['PATNO'] = vol['PATNO'].astype(str)
if 'aseg_eTIV' not in harm.columns:
    harm = harm.merge(vol[['PATNO','Visit','aseg_eTIV']].dropna(), on=['PATNO','Visit'], how='left')

clin_cols = ['MoCA_total','UPDRS3_total','UPDRS2_total','RBDSQ_total','SCOPA_total',
             'GDS_total','ESS_total','UPSIT_total','TrailsB_sec','STAI_state','STAI_trait']
clin_cols = [c for c in clin_cols if c in vol.columns]
harm = harm.drop(columns=[c for c in clin_cols if c in harm.columns], errors='ignore')
harm = harm.merge(vol[['PATNO','Visit']+clin_cols], on=['PATNO','Visit'], how='left')
harm['sex_num'] = (harm['Sex']=='M').astype(float)

mega = safe_read_csv(f'{BASE_DIR}/PPMI_mega_clinical.csv', low_memory=False)
mega.columns = [c.replace('-','_') for c in mega.columns]
mega['PATNO'] = mega['PATNO'].astype(str)
mega['Visit'] = mega['EVENT_ID'].map(EVENT_MAP)

dur = mega[mega['EVENT_ID']=='BL'][['PATNO','disease_duration_at_MRI_yrs']].drop_duplicates('PATNO')
harm = harm.merge(dur, on='PATNO', how='left')

berg = harm[['PATNO','subtype']].dropna().drop_duplicates('PATNO').rename(columns={'subtype':'berg_subtype'})
berg['PATNO'] = berg['PATNO'].astype(str)
patno_berg = berg.set_index('PATNO')['berg_subtype'].to_dict()
harm['subtype'] = harm['PATNO'].map(patno_berg)

_educ_df = safe_read_csv(f'{BASE_DIR}/PPMI_Clinical_13Apr2026/_Subject_Characteristics/Socio-Economics_12Apr2026.csv', low_memory=False)
_educ_df['PATNO'] = _educ_df['PATNO'].astype(str)
_educ_lookup = _educ_df.dropna(subset=['EDUCYRS']).drop_duplicates('PATNO')[['PATNO','EDUCYRS']]
print(f"EDUCYRS available for {len(_educ_lookup)} subjects")
harm = harm.drop(columns=[c for c in harm.columns if c=='EDUCYRS'], errors='ignore')
harm = harm.merge(_educ_lookup, on='PATNO', how='left')

print(f"Loaded: {len(harm)} rows, {harm['PATNO'].nunique()} subjects")
print(f"Groups: {harm['Group'].value_counts().to_dict()}")

amyg_cols  = [c for c in harm.columns if c.startswith('amyg_lh_') or c.startswith('amyg_rh_')]
aparc_cols = [c for c in harm.columns if (c.startswith('lh_aparc_') or c.startswith('rh_aparc_')) and c.endswith('_GrayVol')]
aseg_cols  = [c for c in harm.columns if c.startswith('aseg_') and
              c not in ['aseg_eTIV','aseg_BrainSegVol_to_eTIV','aseg_MaskVol_to_eTIV']]
all_vol = list(set(amyg_cols+aparc_cols+aseg_cols))
all_vol = [c for c in all_vol if harm[c].notna().sum()>100]
qc_cols = [c for c in aparc_cols+aseg_cols if harm[c].notna().sum()>100]

z_df = harm.copy()
for col in all_vol:
    z_df[f'z_{col}'] = z_df[col]
print(f"QC will run on {len(all_vol)} harmonised (unadjusted) volume columns; adjustment follows outlier exclusion")

z_df['visit_rank'] = z_df['Visit'].map(VISIT_ORDER)
fv = z_df.sort_values('visit_rank').groupby('PATNO')['Visit'].first().reset_index().rename(columns={'Visit':'first_visit'})
z_df = z_df.merge(fv, on='PATNO', how='left')
first_mri = z_df[z_df['Visit']==z_df['first_visit']].copy()

qc_z_cols = [f'z_{c}' for c in qc_cols if f'z_{c}' in z_df.columns]
fmri_complete = first_mri.dropna(subset=qc_z_cols)
hc_fmri = fmri_complete[fmri_complete['Group']=='Control']
n_hc = len(hc_fmri)
max_pcs = max(5, n_hc // 5)
qc_scaler_means = hc_fmri[qc_z_cols].mean()
qc_scaler_sds = hc_fmri[qc_z_cols].std().replace(0, 1e-6)
hc_scaled = (hc_fmri[qc_z_cols] - qc_scaler_means) / qc_scaler_sds
all_scaled = (fmri_complete[qc_z_cols] - qc_scaler_means) / qc_scaler_sds
pca = PCA(n_components=max_pcs)
pca.fit(hc_scaled.values)
pca_all = pca.transform(all_scaled.values)
pca_hc = pca.transform(hc_scaled.values)
print(f"HC at first MRI: {n_hc} | PCs: {max_pcs} | Variance explained: {pca.explained_variance_ratio_.sum():.1%}")
mean_hc = pca_hc.mean(axis=0)
cov_inv = np.linalg.inv(np.cov(pca_hc.T))
md_cs = [(fmri_complete.index[i], mahalanobis(row, mean_hc, cov_inv)) for i,row in enumerate(pca_all)]
md_df = pd.DataFrame(md_cs, columns=['idx','MD_cs']).set_index('idx')
threshold_cs = np.sqrt(stats.chi2.ppf(0.999, df=max_pcs))
outlier_patnos_cs = set(z_df.loc[md_df[md_df['MD_cs']>threshold_cs].index, 'PATNO'])
print(f"Cross-sectional outliers: {len(outlier_patnos_cs)} subjects (threshold MD={threshold_cs:.2f})")

_outlier_group_lookup = z_df[z_df['PATNO'].isin(outlier_patnos_cs)][['PATNO','Group']].drop_duplicates('PATNO')
_outlier_group_counts = _outlier_group_lookup['Group'].value_counts()
print(f"  By group: {_outlier_group_counts.to_dict()}")
if 'subtype' in z_df.columns:
    _outlier_subtype_lookup = z_df[(z_df['PATNO'].isin(outlier_patnos_cs)) & (z_df['Group']=='PD')][['PATNO','subtype']].drop_duplicates('PATNO')
    _outlier_subtype_counts = _outlier_subtype_lookup['subtype'].value_counts(dropna=False)
    print(f"  PD outliers by subtype: {_outlier_subtype_counts.to_dict()}")

amyg_z_cols = [f'z_{c}' for c in amyg_cols if f'z_{c}' in z_df.columns]
long_df = z_df[z_df['Visit'].isin(VISITS)].copy()
outlier_slope_patnos = set()
for col in amyg_z_cols:
    for grp_df in [long_df[long_df['Group']=='Control'], long_df[long_df['Group']=='PD']]:
        slopes = {}
        for patno, g in grp_df.groupby('PATNO'):
            tmp = g[[col,'years_from_first_mri']].dropna()
            if len(tmp)<2: continue
            s,_,_,_,_ = stats.linregress(tmp['years_from_first_mri'], tmp[col])
            slopes[patno] = s
        if len(slopes)<10: continue
        slopes = pd.Series(slopes)
        med, mad = slopes.median(), (slopes-slopes.median()).abs().median()*1.4826
        if mad<1e-6: continue
        outlier_slope_patnos.update(slopes[(slopes-med).abs()/mad>4].index.tolist())
print(f"Longitudinal trajectory outliers: {len(outlier_slope_patnos)} subjects")

_outlier_slope_group_lookup = z_df[z_df['PATNO'].isin(outlier_slope_patnos)][['PATNO','Group']].drop_duplicates('PATNO')
_outlier_slope_group_counts = _outlier_slope_group_lookup['Group'].value_counts()
print(f"  By group: {_outlier_slope_group_counts.to_dict()}")
if 'subtype' in z_df.columns:
    _outlier_slope_subtype_lookup = z_df[(z_df['PATNO'].isin(outlier_slope_patnos)) & (z_df['Group']=='PD')][['PATNO','subtype']].drop_duplicates('PATNO')
    _outlier_slope_subtype_counts = _outlier_slope_subtype_lookup['subtype'].value_counts(dropna=False)
    print(f"  PD outliers by subtype: {_outlier_slope_subtype_counts.to_dict()}")

all_outliers = outlier_patnos_cs | outlier_slope_patnos
df = z_df[~z_df['PATNO'].isin(all_outliers)].copy()
_adj_cols = ['Age','sex_num','aseg_eTIV']
for col in all_vol:
    _zc = f'z_{col}'
    df[_zc] = np.nan
    _tmp = df[_adj_cols+[col]].dropna()
    if len(_tmp) < 20: continue
    _X = sm.add_constant(_tmp[_adj_cols])
    _m = sm.OLS(_tmp[col], _X).fit()
    df.loc[_tmp.index, _zc] = (_tmp[col] - _X.dot(_m.params)) + _tmp[col].mean()
print(f"Adjusted {len(all_vol)} volume columns AFTER outlier exclusion (eTIV, age, sex; whole retained sample; grand mean added back, mm3)")
print(f"Total excluded: {len(all_outliers)} subjects | Remaining: {df['PATNO'].nunique()} subjects, {len(df)} rows")

_total_outlier_group_lookup = z_df[z_df['PATNO'].isin(all_outliers)][['PATNO','Group']].drop_duplicates('PATNO')
_total_outlier_group_counts = _total_outlier_group_lookup['Group'].value_counts()
print(f"  Total excluded, by group: {_total_outlier_group_counts.to_dict()}")
_qc_group_rows = []
for grp_name, grp_count in _total_outlier_group_counts.items():
    _qc_group_rows.append({'exclusion_type':'total','group':grp_name,'n_excluded':grp_count})
for grp_name, grp_count in _outlier_group_counts.items():
    _qc_group_rows.append({'exclusion_type':'cross_sectional','group':grp_name,'n_excluded':grp_count})
for grp_name, grp_count in _outlier_slope_group_counts.items():
    _qc_group_rows.append({'exclusion_type':'longitudinal','group':grp_name,'n_excluded':grp_count})
pd.DataFrame(_qc_group_rows).to_csv(f'{OUT_DIR}/qc_outlier_group_breakdown.csv', index=False)
print(f"  Saved: {OUT_DIR}/qc_outlier_group_breakdown.csv")

df.to_csv(f'{OUT_DIR}/harmonised_adjustedvol_qc.csv', index=False)
