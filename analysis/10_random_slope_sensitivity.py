"""Random-slope sensitivity models."""

def _fit_tiered_9j(data, col, group_col, time_col='years_from_first_mri', is_pd_hc=False):
    if is_pd_hc:
        sub = data[data['Visit'].isin(VISITS)].dropna(subset=[col, time_col, group_col]).copy()
    else:
        sub = data[(data['Group']=='PD') & (data['Visit'].isin(VISITS))].dropna(subset=[col, time_col, group_col]).copy()
    multi = sub.groupby('PATNO').filter(lambda x: len(x) > 1)
    if multi['PATNO'].nunique() < 15: return None, 'insufficient_data'
    col_s = col.replace('-','_'); multi_c = multi.copy(); multi_c[col_s] = multi_c[col]
    formula = f'{col_s} ~ {time_col} * {group_col}'

    try:
        mdl = smf.mixedlm(formula, data=multi_c, groups=multi_c['PATNO'], re_formula=f'~{time_col}').fit(reml=True)
        if mdl.converged:
            return mdl, 'random_slope_full'
    except Exception:
        pass

    try:
        vc = {'time_re': f'0 + {time_col}'}
        mdl = smf.mixedlm(formula, data=multi_c, groups=multi_c['PATNO'], re_formula='1', vc_formula=vc).fit(reml=True)
        if mdl.converged:
            return mdl, 'random_slope_zerocorr'
    except Exception:
        pass

    try:
        mdl = smf.mixedlm(formula, data=multi_c, groups=multi_c['PATNO']).fit(reml=True)
        if mdl.converged:
            return mdl, 'random_intercept_only'
    except Exception:
        pass

    return None, 'failed_all_tiers'

def _report_tiered_interaction_9j(mdl, tier, group_col, time_col='years_from_first_mri'):
    if mdl is None:
        return None, None
    int_key = [k for k in mdl.params.index if group_col in k and time_col in k and ':' in k]
    if not int_key:
        return None, None
    return mdl.params[int_key[0]], mdl.pvalues[int_key[0]]

print("\nPD vs HC: tiered random-slope sensitivity check")
_tier_summary_9j = []
_df_9j = df.copy()
_df_9j['is_PD'] = (_df_9j['Group']=='PD').astype(float)
for col, wlabel in [(whole_col, whole_label)] + SECTORS:
    mdl, tier = _fit_tiered_9j(_df_9j, col, 'is_PD', is_pd_hc=True)
    beta, p = _report_tiered_interaction_9j(mdl, tier, 'is_PD')
    if beta is not None:
        print(f"  {wlabel:<16}: tier={tier:<22} interaction={beta:+.3f} p={p:.4f}{sig(p)}")
    else:
        print(f"  {wlabel:<16}: tier={tier} -- no interaction term found")
    _tier_summary_9j.append({'analysis':'PD_vs_HC','region':wlabel,'tier':tier,'interaction':beta,'p':p})

print("\nEarly/non-early VH: tiered random-slope sensitivity check")
_vh_valid_9j = vh_early_groups.dropna(subset=['early']).rename(columns={'early':'early_VH_9j'})
_vh_long_9j = df.merge(_vh_valid_9j[['PATNO','early_VH_9j']], on='PATNO', how='inner')
for col, wlabel in [(whole_col, whole_label)] + SECTORS:
    mdl, tier = _fit_tiered_9j(_vh_long_9j, col, 'early_VH_9j')
    beta, p = _report_tiered_interaction_9j(mdl, tier, 'early_VH_9j')
    if beta is not None:
        print(f"  {wlabel:<16}: tier={tier:<22} interaction={beta:+.3f} p={p:.4f}{sig(p)}")
    else:
        print(f"  {wlabel:<16}: tier={tier} -- no interaction term found")
    _tier_summary_9j.append({'analysis':'Early_Late_VH','region':wlabel,'tier':tier,'interaction':beta,'p':p})

print("\nEarly/non-early CI: tiered random-slope sensitivity check")
_ci_valid_9j = ci_early_groups.dropna(subset=['early']).rename(columns={'early':'early_CI_9j'})
_ci_long_9j = df.merge(_ci_valid_9j[['PATNO','early_CI_9j']], on='PATNO', how='inner')
for col, wlabel in [(whole_col, whole_label)] + SECTORS:
    mdl, tier = _fit_tiered_9j(_ci_long_9j, col, 'early_CI_9j')
    beta, p = _report_tiered_interaction_9j(mdl, tier, 'early_CI_9j')
    if beta is not None:
        print(f"  {wlabel:<16}: tier={tier:<22} interaction={beta:+.3f} p={p:.4f}{sig(p)}")
    else:
        print(f"  {wlabel:<16}: tier={tier} -- no interaction term found")
    _tier_summary_9j.append({'analysis':'Early_Late_CI','region':wlabel,'tier':tier,'interaction':beta,'p':p})

pd.DataFrame(_tier_summary_9j).to_csv(f'{OUT_DIR}/Step9J_random_slope_sensitivity.csv', index=False)
print(f"\nSaved: {OUT_DIR}/Step9J_random_slope_sensitivity.csv")
