"""PD versus HC, cross-sectional and longitudinal."""

def welch_t(hc_v, pd_v):
    t, p = stats.ttest_ind(hc_v, pd_v, equal_var=False)
    return t, p

whole_col, whole_label = WHOLE
hc_v = get_grp(bl_all,'HC')[whole_col].dropna().values
pd_v = get_grp(bl_all,'All_PD')[whole_col].dropna().values
t_whole, p_whole = welch_t(hc_v, pd_v)
print(f"\nWhole amygdala (standalone t-test)")
print(f"  {whole_label:<16}: HC={hc_v.mean():.1f} | PD={pd_v.mean():.1f} mm3 | t={t_whole:.3f} p={p_whole:.4f}{sig(p_whole)} "
      f"n_HC={len(hc_v)} n_PD={len(pd_v)}")
pd.DataFrame([{'region':whole_label,'col':whole_col,'mean_HC_mm3':hc_v.mean(),'mean_PD_mm3':pd_v.mean(),
              't':t_whole,'p':p_whole,'n_HC':len(hc_v),'n_PD':len(pd_v)}]).to_csv(f'{OUT_DIR}/cs4_whole_amygdala_PD_vs_HC.csv', index=False)

print(f"\nSectors sCLR/CMA/BLA (t-tests, FDR family of 3)")
sec_rows = []; sec_ps = []; sec_meta = []
for col, label in SECTORS:
    hc_v = get_grp(bl_all,'HC')[col].dropna().values
    pd_v = get_grp(bl_all,'All_PD')[col].dropna().values
    t, p = welch_t(hc_v, pd_v)
    sec_ps.append(p); sec_meta.append((col,label,hc_v,pd_v,t))
sec_p_fdr = multipletests(sec_ps, method='fdr_bh')[1] if sec_ps else []
for (col,label,hc_v,pd_v,t), p, p_fdr in zip(sec_meta, sec_ps, sec_p_fdr):
    print(f"  {label:<8}: HC={hc_v.mean():.1f} | PD={pd_v.mean():.1f} mm3 | t={t:.3f} p={p:.4f} p_fdr={p_fdr:.4f}{sig(p_fdr)}")
    sec_rows.append({'region':label,'col':col,'mean_HC_mm3':hc_v.mean(),'mean_PD_mm3':pd_v.mean(),
                     't':t,'p':p,'p_fdr':p_fdr,'n_HC':len(hc_v),'n_PD':len(pd_v)})
pd.DataFrame(sec_rows).to_csv(f'{OUT_DIR}/cs4_sectors_PD_vs_HC.csv', index=False)

def lme_group_fixed_slope(data, col, group_val, time_col='years_from_first_mri'):
    sub = get_grp(data[data['Visit'].isin(VISITS)], group_val)
    sub = sub.dropna(subset=[col, time_col]).copy()
    multi = sub.groupby('PATNO').filter(lambda x: len(x) > 1)
    if multi['PATNO'].nunique() < 8: return None
    col_s = col.replace('-','_')
    multi = multi.copy(); multi[col_s] = multi[col]
    try:
        mdl = smf.mixedlm(f'{col_s} ~ {time_col}', data=multi, groups=multi['PATNO']).fit(reml=True)
        return {'slope':mdl.params[time_col], 'p':mdl.pvalues[time_col], 'n':multi['PATNO'].nunique()}
    except: return None

def lme_group_time_interaction(data, col, time_col='years_from_first_mri'):
    sub = data[data['Visit'].isin(VISITS)].copy()
    sub = sub[sub['Group'].isin(['Control','PD'])].copy()
    sub = sub.dropna(subset=[col, time_col]).copy()
    multi = sub.groupby('PATNO').filter(lambda x: len(x) > 1)
    if multi['PATNO'].nunique() < 15: return None
    col_s = col.replace('-','_')
    multi = multi.copy(); multi[col_s] = multi[col]
    multi['is_PD'] = (multi['Group']=='PD').astype(float)
    try:
        mdl = smf.mixedlm(f'{col_s} ~ {time_col} * is_PD', data=multi, groups=multi['PATNO']).fit(reml=True)
        int_key = [k for k in mdl.params.index if 'is_PD' in k and time_col in k and ':' in k]
        int_term = int_key[0] if int_key else None
        if int_term is None: return None
        return {'n':multi['PATNO'].nunique(), 'interaction_beta':mdl.params[int_term],
                'p_interaction':mdl.pvalues[int_term], 'converged':mdl.converged}
    except: return None

print("\nFixed per-group slopes: HC vs All_PD")
fixed_rows = []
for col, label in [WHOLE] + SECTORS:
    for grp in ['HC','All_PD']:
        res = lme_group_fixed_slope(df, col, grp)
        if res is None: continue
        fixed_rows.append({'region':label,'col':col,'group':grp, **res})
_fx = pd.DataFrame(fixed_rows)
if len(_fx):
    _fx['p_fdr'] = np.nan
    for _g5 in _fx['group'].unique():
        _m5 = (_fx['group'] == _g5) & (_fx['col'] != WHOLE[0])
        if _m5.sum() > 0:
            _fx.loc[_m5, 'p_fdr'] = multipletests(_fx.loc[_m5, 'p'].values, method='fdr_bh')[1]
    for _r5 in _fx.itertuples():
        if _r5.col == WHOLE[0]:
            print(f"  {_r5.region:<16} | {_r5.group:<7}: slope={_r5.slope:+.3f} mm3/yr p={_r5.p:.4f}{sig(_r5.p)} n={_r5.n}")
        else:
            print(f"  {_r5.region:<16} | {_r5.group:<7}: slope={_r5.slope:+.3f} mm3/yr p={_r5.p:.4f} p_fdr={_r5.p_fdr:.4f}{sig(_r5.p_fdr)} n={_r5.n}")
_fx.to_csv(f'{OUT_DIR}/long5a_fixed_slopes_PD_vs_HC.csv', index=False)

print("\nGroup x time interaction: whole amygdala standalone")
int_whole = lme_group_time_interaction(df, whole_col)
int_whole_rows = []
if int_whole is not None:
    print(f"  {whole_label:<16}: interaction={int_whole['interaction_beta']:+.3f} p={int_whole['p_interaction']:.4f}{sig(int_whole['p_interaction'])} n={int_whole['n']}")
    int_whole_rows.append({'region':whole_label,'col':whole_col, **int_whole})
pd.DataFrame(int_whole_rows).to_csv(f'{OUT_DIR}/long5b_whole_amygdala_interaction.csv', index=False)

print("\nGroup x time interaction: sectors sCLR/CMA/BLA (FDR family of 3)")
int_sec_rows = []; int_sec_ps = []; int_sec_meta = []
for col, label in SECTORS:
    res = lme_group_time_interaction(df, col)
    if res is None: continue
    int_sec_ps.append(res['p_interaction']); int_sec_meta.append((col,label,res))
int_sec_p_fdr = multipletests(int_sec_ps, method='fdr_bh')[1] if int_sec_ps else []
for (col,label,res), p_fdr in zip(int_sec_meta, int_sec_p_fdr):
    print(f"  {label:<8}: interaction={res['interaction_beta']:+.3f} p={res['p_interaction']:.4f} p_fdr={p_fdr:.4f}{sig(p_fdr)} n={res['n']}")
    int_sec_rows.append({'region':label,'col':col,'p_fdr':p_fdr, **res})
pd.DataFrame(int_sec_rows).to_csv(f'{OUT_DIR}/long5b_sectors_interaction.csv', index=False)

any_sector_interaction_sig = (pd.DataFrame(int_sec_rows)['p_fdr'] < 0.05).any() if int_sec_rows else False
print(f"At least one sector Group x time interaction FDR-significant: {any_sector_interaction_sig}")
nuc_interaction_rows = []
if any_sector_interaction_sig:
    print("\nPost-hoc: nucleus-level Group x time interaction (own FDR family of 9)")
    nuc_int_ps = []; nuc_int_meta = []
    for col, label in NUCLEI:
        res = lme_group_time_interaction(df, col)
        if res is None: continue
        nuc_int_ps.append(res['p_interaction']); nuc_int_meta.append((col,label,res))
    nuc_int_p_fdr = multipletests(nuc_int_ps, method='fdr_bh')[1] if nuc_int_ps else []
    for (col,label,res), p_fdr in zip(nuc_int_meta, nuc_int_p_fdr):
        print(f"  {label:<28}: interaction={res['interaction_beta']:+.3f} p={res['p_interaction']:.4f} "
              f"p_fdr={p_fdr:.4f}{sig(p_fdr)} n={res['n']}")
        nuc_interaction_rows.append({'region':label,'col':col,'p_fdr':p_fdr, **res})
else:
    print("\n  No sector cleared FDR significance -- nucleus-level post-hoc skipped.")
pd.DataFrame(nuc_interaction_rows).to_csv(f'{OUT_DIR}/long5b_nuclei_interaction_posthoc.csv', index=False)
