"""Early versus non-early VH and CI, cross-sectional, longitudinal, MoCA-adjusted and demographic comparisons."""

_ALL_NUCLEI_9D = ['Lateral_nucleus','Basal_nucleus','Accessory_Basal_nucleus','Paralaminar_nucleus',
                  'Central_nucleus','Medial_nucleus','Cortical_nucleus','Anterior_amygdaloid_area_AAA','Corticoamygdaloid_transitio']
for _n in _ALL_NUCLEI_9D:
    _lz, _rz = f'z_amyg_lh_{_n}', f'z_amyg_rh_{_n}'
    if _lz in df.columns and _rz in df.columns and f'z_bi_{_n}' not in df.columns:
        df[f'z_bi_{_n}'] = (df[_lz] + df[_rz]) / 2
_NUCLEI_COLS_9D = [(f'z_bi_{_n}', _n.replace('_',' ')) for _n in _ALL_NUCLEI_9D if f'z_bi_{_n}' in df.columns]

def _welch_t_d_9d(pos, neg):
    t, p = stats.ttest_ind(pos, neg, equal_var=False)
    d = (pos.mean()-neg.mean())/np.sqrt((pos.std()**2+neg.std()**2)/2+1e-9)
    return t, p, d

def _run_cross_sectional_9d(groups_df, group_col, label):
    valid = groups_df.dropna(subset=['early']).rename(columns={'early': group_col})
    bl_data = bl_pd.merge(valid[['PATNO', group_col]], on='PATNO', how='inner')
    print(f"\n{label} cross-sectional: n={len(bl_data)}, early={int(bl_data[group_col].sum())}, not_early={int((bl_data[group_col]==0).sum())}")

    pos = bl_data[bl_data[group_col]==1][whole_col].dropna()
    neg = bl_data[bl_data[group_col]==0][whole_col].dropna()
    if len(pos)>=3 and len(neg)>=3:
        t,p,d = _welch_t_d_9d(pos, neg)
        print(f"  {whole_label:<16} (standalone): d={d:+.3f} p={p:.4f}{sig(p)} n_early={len(pos)} n_not={len(neg)}")

    sec_results = []
    for col, wlabel in SECTORS:
        pos = bl_data[bl_data[group_col]==1][col].dropna()
        neg = bl_data[bl_data[group_col]==0][col].dropna()
        if len(pos)<3 or len(neg)<3: continue
        t,p,d = _welch_t_d_9d(pos, neg)
        sec_results.append({'region':wlabel,'col':col,'d':d,'p':p,'n_pos':len(pos),'n_neg':len(neg)})
    if sec_results:
        p_fdr = multipletests([r['p'] for r in sec_results], method='fdr_bh')[1]
        for r, pf in zip(sec_results, p_fdr):
            r['p_fdr'] = pf
            print(f"  {r['region']:<8}: d={r['d']:+.3f} p={r['p']:.4f} p_fdr={pf:.4f}{sig(pf)} n_early={r['n_pos']} n_not={r['n_neg']}")

    any_sig = any(r.get('p_fdr',1) < 0.05 for r in sec_results)
    print(f"  At least one sector significant: {any_sig}")

    nuc_results = []
    if any_sig:
        print(f"\n  {label} nucleus post-hoc (own FDR family of 9)")
        for col, nlabel in _NUCLEI_COLS_9D:
            pos = bl_data[bl_data[group_col]==1][col].dropna()
            neg = bl_data[bl_data[group_col]==0][col].dropna()
            if len(pos)<3 or len(neg)<3: continue
            t,p,d = _welch_t_d_9d(pos, neg)
            nuc_results.append({'region':nlabel,'col':col,'d':d,'p':p,'n_pos':len(pos),'n_neg':len(neg)})
        if nuc_results:
            p_fdr_n = multipletests([r['p'] for r in nuc_results], method='fdr_bh')[1]
            for r, pf in zip(nuc_results, p_fdr_n):
                r['p_fdr'] = pf
                print(f"    {r['region']:<28}: d={r['d']:+.3f} p={r['p']:.4f} p_fdr={pf:.4f}{sig(pf)}")
    pd.DataFrame(sec_results).to_csv(f'{OUT_DIR}/{label}_9D_sectors.csv', index=False)
    pd.DataFrame(nuc_results).to_csv(f'{OUT_DIR}/{label}_9D_nuclei.csv', index=False)
    return bl_data, sec_results, nuc_results

vh_bl_data_9d, vh_sec_results_9d, vh_nuc_results_9d = _run_cross_sectional_9d(vh_early_groups, 'early_VH', 'VH')
ci_bl_data_9d, ci_sec_results_9d, ci_nuc_results_9d = _run_cross_sectional_9d(ci_early_groups, 'early_CI', 'CI')

def _within_fdr_9e(rows):
    from statsmodels.stats.multitest import multipletests as _mtw
    d_ = pd.DataFrame(rows)
    if d_.empty: return d_
    d_['p_fdr'] = np.nan
    for g_ in d_['group'].unique():
        m_ = (d_['group'] == g_) & (~d_['whole'])
        if m_.sum() > 0:
            d_.loc[m_, 'p_fdr'] = _mtw(d_.loc[m_, 'p'].values, method='fdr_bh')[1]
    for r_ in d_.itertuples():
        tier_ = f" tier={r_.tier}" if hasattr(r_, 'tier') else ''
        if r_.whole:
            print(f"  {r_.region:<16} | {r_.group:<10}: slope={r_.slope:+.3f} p={r_.p:.4f}{sig(r_.p)} n={r_.n}{tier_}")
        else:
            print(f"  {r_.region:<16} | {r_.group:<10}: slope={r_.slope:+.3f} p={r_.p:.4f} p_fdr={r_.p_fdr:.4f}{sig(r_.p_fdr)} n={r_.n}{tier_}")
    return d_

def _fit_within_group_slope_9e(data, col, group_col, group_val, time_col='years_from_first_mri'):
    sub = data[(data['Group']=='PD') & (data['Visit'].isin(VISITS)) & (data[group_col]==group_val)].dropna(subset=[col, time_col]).copy()
    multi = sub.groupby('PATNO').filter(lambda x: len(x) > 1)
    if multi['PATNO'].nunique() < 10: return None
    col_s = col.replace('-','_'); multi_c = multi.copy(); multi_c[col_s] = multi_c[col]
    try:
        return smf.mixedlm(f'{col_s} ~ {time_col}', data=multi_c, groups=multi_c['PATNO']).fit(reml=True), multi['PATNO'].nunique()
    except: return None

def _fit_time_model_9e(data, col, group_col, time_col='years_from_first_mri'):
    sub = data[(data['Group']=='PD') & (data['Visit'].isin(VISITS))].dropna(subset=[col, time_col, group_col]).copy()
    multi = sub.groupby('PATNO').filter(lambda x: len(x) > 1)
    if multi['PATNO'].nunique() < 15: return None
    gcounts = multi.drop_duplicates('PATNO')[group_col].value_counts()
    if (gcounts < 10).any(): return None
    col_s = col.replace('-','_'); multi_c = multi.copy(); multi_c[col_s] = multi_c[col]
    try:
        return smf.mixedlm(f'{col_s} ~ {time_col} * {group_col}', data=multi_c, groups=multi_c['PATNO']).fit(reml=True)
    except: return None

def _run_longitudinal_9e(groups_df, group_col, label):
    valid = groups_df.dropna(subset=['early']).rename(columns={'early': group_col})
    long_data = df.merge(valid[['PATNO', group_col]], on='PATNO', how='inner')

    print(f"\n{label}: within-group time effects (fixed per-group slopes)")
    _wg_rows = []
    for col, wlabel in [(whole_col, whole_label)] + SECTORS:
        for gval, glabel in [(1,'early'), (0,'not_early')]:
            res = _fit_within_group_slope_9e(long_data, col, group_col, gval)
            if res is None:
                print(f"  {wlabel:<16} | {glabel:<10}: insufficient data"); continue
            mdl, n = res
            time_key = [k for k in mdl.params.index if k != 'Intercept'][0]
            slope = mdl.params[time_key]; p = mdl.pvalues[time_key]
            _wg_rows.append({'region':wlabel,'group':glabel,'slope':slope,'p':p,'n':n,'whole':(col == whole_col)})
    _wg_df = _within_fdr_9e(_wg_rows)
    if len(_wg_df): _wg_df.drop(columns='whole').to_csv(f'{OUT_DIR}/{label}_9E_within_group_slopes.csv', index=False)

    print(f"\n{label}: group x time interaction (FDR-corrected)")
    _fitted_mdls_9e = {}
    for col, wlabel in [(whole_col, whole_label)] + SECTORS:
        mdl = _fit_time_model_9e(long_data, col, group_col)
        _fitted_mdls_9e[col] = mdl

    _mdl_w9e = _fitted_mdls_9e[whole_col]
    if _mdl_w9e is None:
        print(f"  {whole_label}: insufficient data")
    else:
        _int_key_w9e = [k for k in _mdl_w9e.params.index if group_col in k and 'years_from_first_mri' in k and ':' in k]
        if _int_key_w9e:
            print(f"  {whole_label:<16} (standalone): interaction={_mdl_w9e.params[_int_key_w9e[0]]:+.3f} p={_mdl_w9e.pvalues[_int_key_w9e[0]]:.4f}{sig(_mdl_w9e.pvalues[_int_key_w9e[0]])}")

    _sec_int_ps_9e = []; _sec_int_meta_9e = []
    for col, wlabel in SECTORS:
        mdl = _fitted_mdls_9e[col]
        if mdl is None:
            print(f"  {wlabel}: insufficient data"); continue
        int_key = [k for k in mdl.params.index if group_col in k and 'years_from_first_mri' in k and ':' in k]
        if int_key:
            _sec_int_ps_9e.append(mdl.pvalues[int_key[0]])
            _sec_int_meta_9e.append((wlabel, mdl.params[int_key[0]], mdl.pvalues[int_key[0]]))
    if _sec_int_ps_9e:
        _sec_p_fdr_9e = multipletests(_sec_int_ps_9e, method='fdr_bh')[1]
        for (wlabel, beta, p), pf in zip(_sec_int_meta_9e, _sec_p_fdr_9e):
            print(f"  {wlabel:<16}: interaction={beta:+.3f} p={p:.4f} p_fdr={pf:.4f}{sig(pf)}")

    print(f"\n{label}: EMM by timepoint (years 0-4), FDR-corrected contrasts")
    _emm_cache = {}
    for col, wlabel in [(whole_col, whole_label)] + SECTORS:
        mdl = _fit_time_model_9e(long_data, col, group_col)
        if mdl is None:
            print(f"  {wlabel}: insufficient data"); continue
        time_col = 'years_from_first_mri'
        contrast_ps = []
        rows_this_region = []
        for t in [0,1,2,3,4]:
            p_params = mdl.params; cov = mdl.cov_params()
            b0=p_params.get('Intercept',0.0); b_t=p_params.get(time_col,0.0); b_g=p_params.get(group_col,0.0)
            int_key2 = [k for k in p_params.index if group_col in k and time_col in k and ':' in k]
            b_int = p_params.get(int_key2[0],0.0) if int_key2 else 0.0
            def _emm(is_pos):
                x_t,x_g,x_int = t, float(is_pos), t*float(is_pos)
                emm = b0+b_t*x_t+b_g*x_g+b_int*x_int
                names = ['Intercept',time_col,group_col]+([int_key2[0]] if int_key2 else [])
                a = np.array([1.0,x_t,x_g]+([x_int] if int_key2 else []))
                try:
                    C = cov.loc[names,names].values
                    se = np.sqrt(max(a@C@a,0))
                except: se=np.nan
                return emm, se
            e0, se0 = _emm(0); e1, se1 = _emm(1)
            names_c = [group_col]+([int_key2[0]] if int_key2 else [])
            a_c = np.array([1.0]+([t] if int_key2 else []))
            diff = b_g + t*b_int
            try:
                C_c = cov.loc[names_c,names_c].values
                se_d = np.sqrt(max(a_c@C_c@a_c,0))
                z = diff/se_d if se_d>0 else np.nan
                pv = 2*(1-stats.norm.cdf(abs(z))) if not np.isnan(z) else np.nan
            except: se_d,pv = np.nan,np.nan
            contrast_ps.append(pv)
            rows_this_region.append((t,e0,se0,e1,se1,diff,se_d,pv))
        p_fdr_arr = multipletests([p if not np.isnan(p) else 1.0 for p in contrast_ps], method='fdr_bh')[1]
        print(f"  {wlabel}:")
        for (t,e0,se0,e1,se1,diff,se_d,pv), pf in zip(rows_this_region, p_fdr_arr):
            print(f"    Year {t}: not_early={e0:.2f}(SE={se0:.2f}) early={e1:.2f}(SE={se1:.2f}) diff={diff:+.2f} p={pv:.4f} p_fdr={pf:.4f}{sig(pf)}")
        _emm_cache[col] = {
            'years': [r[0] for r in rows_this_region],
            'e0s': [r[1] for r in rows_this_region], 'se0s': [r[2] for r in rows_this_region],
            'e1s': [r[3] for r in rows_this_region], 'se1s': [r[4] for r in rows_this_region],
            'p_fdr': list(p_fdr_arr),
            'p_uncorrected': [r[7] for r in rows_this_region],
        }
    return long_data, _emm_cache

vh_long_data_9e, vh_emm_cache_9e = _run_longitudinal_9e(vh_early_groups, 'early_VH', 'VH')
ci_long_data_9e, ci_emm_cache_9e = _run_longitudinal_9e(ci_early_groups, 'early_CI', 'CI')

_moca_col_9f = [c for c in mega.columns if c.upper() == 'MOCA_MCATOT']
if _moca_col_9f:
    _moca_visits_9f = mega[mega['PATNO'].notna() & mega[_moca_col_9f[0]].notna()][['PATNO','Visit',_moca_col_9f[0]]].copy()
    _moca_visits_9f['visit_rank'] = _moca_visits_9f['Visit'].map(VISIT_ORDER)
    _moca_visits_9f = _moca_visits_9f.dropna(subset=['visit_rank'])
    _moca_lookup_9f = _moca_visits_9f.sort_values('visit_rank').drop_duplicates('PATNO', keep='first')[['PATNO', _moca_col_9f[0]]]
    print(f"MoCA at earliest available visit: {len(_moca_lookup_9f)} participants")
else:
    _moca_lookup_9f = None
_moca_col_name_9f = _moca_col_9f[0] if _moca_col_9f else None

def _run_moca_sensitivity_9f(bl_data, group_col, label):
    if _moca_col_name_9f is None:
        print("MoCA column not found -- skipping"); return
    data = bl_data.merge(_moca_lookup_9f, on='PATNO', how='left') if _moca_col_name_9f not in bl_data.columns else bl_data
    print(f"\n{label}: MoCA-adjusted sensitivity")
    for col, wlabel in [('z_sCLR','sCLR'), ('z_CMA','CMA')]:
        for model_label, covs in [('Unadjusted', []), ('+ MoCA', [_moca_col_name_9f])]:
            cols_needed = [group_col, col] + covs
            tmp = data.dropna(subset=cols_needed).copy()
            if len(tmp) < 20: continue
            X = sm.add_constant(tmp[[col] + covs])
            try:
                mdl = sm.Logit(tmp[group_col], X).fit(disp=0)
                OR = np.exp(mdl.params[col]); pv = mdl.pvalues[col]
                ci_lo, ci_hi = np.exp(mdl.conf_int().loc[col])
                OR_vol = 1/OR; ci_lo_vol, ci_hi_vol = 1/ci_hi, 1/ci_lo
                print(f"  {wlabel} {model_label:<12}: OR(higher volume)={OR_vol:.3f} [{ci_lo_vol:.3f}-{ci_hi_vol:.3f}] p={pv:.4f}{sig(pv)} n={len(tmp)}")
            except Exception as e:
                print(f"  {wlabel} {model_label:<12}: FAILED ({e})")

_run_moca_sensitivity_9f(vh_bl_data_9d, 'early_VH', 'VH')
_run_moca_sensitivity_9f(ci_bl_data_9d, 'early_CI', 'CI')

_EXTRA_VARS_9H = {
    'Education (years)': 'EDUCYRS',
    'TMT-B (seconds)': 'TrailsB_sec',
    'MDS-UPDRS III': 'UPDRS3_total',
    'MDS-UPDRS II': 'UPDRS2_total',
    'RBDSQ': 'RBDSQ_total',
    'SCOPA-AUT': 'SCOPA_total',
    'GDS': 'GDS_total',
    'ESS': 'ESS_total',
    'UPSIT': 'UPSIT_total',
}

def _run_demographics_9h(bl_data, group_col, label):
    print(f"\n{label}: demographics, early vs not-early")

    _cog_test_cols_9h = ['HVLT_DVT_TOTAL_RECALL','HVLT_DVT_RECOG_DISC_INDEX','SEMFLU_DVS_SFTANIM',
                          'SDMT_DVT_SDM','LNS_DVS_LNS']
    _cog_present_9h = [c for c in _cog_test_cols_9h if c in mega.columns]
    _all_merge_cols_9h = list(set(_cog_present_9h))

    if _all_merge_cols_9h:
        _v9h = mega[mega['PATNO'].notna()][['PATNO','Visit']+_all_merge_cols_9h].copy()
        _v9h['visit_rank'] = _v9h['Visit'].map(VISIT_ORDER)
        _v9h = _v9h.dropna(subset=['visit_rank'])
        _lookup9h = _v9h.sort_values('visit_rank').drop_duplicates('PATNO', keep='first')[['PATNO']+_all_merge_cols_9h]
        bl_data = bl_data.merge(_lookup9h, on='PATNO', how='left', suffixes=('','_dup9h'))

    if 'EDUCYRS' not in bl_data.columns:
        bl_data = bl_data.merge(_educ_lookup, on='PATNO', how='left')

    if 'DVS_JLO_MSSA' not in bl_data.columns:
        try:
            _jlo_9h = safe_read_csv(f'{BASE_DIR}/PPMI_Clinical_13Apr2026/Non-motor_Assessments/Benton_Judgement_of_Line_Orientation_12Apr2026.csv', low_memory=False)
            _jlo_9h.columns = [c.replace('-','_') for c in _jlo_9h.columns]
            _jlo_9h['PATNO'] = _jlo_9h['PATNO'].astype(str)
            _jlo_9h['Visit'] = _jlo_9h['EVENT_ID'].map(EVENT_MAP)
            _jlo_9h['visit_rank'] = _jlo_9h['Visit'].map(VISIT_ORDER)
            _jlo_9h = _jlo_9h.dropna(subset=['visit_rank','DVS_JLO_MSSA'])
            _jlo_lookup_9h = _jlo_9h.sort_values('visit_rank').drop_duplicates('PATNO', keep='first')[['PATNO','DVS_JLO_MSSA']]
            bl_data = bl_data.merge(_jlo_lookup_9h, on='PATNO', how='left')
        except Exception as e:
            print(f"  Could not load Benton JLO: {e}")

    moca_c = [c for c in bl_data.columns if 'MOCA' in c.upper() and 'TOTAL' in c.upper()]

    table_spec = [('Age (years)', 'Age'), ('Education (years)', 'EDUCYRS'),
                  ('Disease duration (years)', 'disease_duration_at_MRI_yrs')]
    table_spec += [('MoCA total', moca_c[0] if moca_c else None)]
    table_spec += [('TMT-B (seconds)', _EXTRA_VARS_9H['TMT-B (seconds)'])]
    table_spec += [('MDS-UPDRS III', _EXTRA_VARS_9H['MDS-UPDRS III']), ('MDS-UPDRS II', _EXTRA_VARS_9H['MDS-UPDRS II'])]
    table_spec += [('RBDSQ', _EXTRA_VARS_9H['RBDSQ']), ('SCOPA-AUT', _EXTRA_VARS_9H['SCOPA-AUT']),
                    ('GDS', _EXTRA_VARS_9H['GDS']), ('ESS', _EXTRA_VARS_9H['ESS']), ('UPSIT', _EXTRA_VARS_9H['UPSIT'])]
    table_spec += [('HVLT total recall', 'HVLT_DVT_TOTAL_RECALL'), ('HVLT recognition discrimination', 'HVLT_DVT_RECOG_DISC_INDEX'),
                    ('Semantic fluency', 'SEMFLU_DVS_SFTANIM'), ('SDMT', 'SDMT_DVT_SDM'), ('LNS', 'LNS_DVS_LNS'),
                    ('Benton JLO', 'DVS_JLO_MSSA')]

    table_rows = []
    print(f"\n{'Variable':<32}{'Early':<22}{'Non-early':<22}{'p':<10}")
    print("-"*86)
    for var_label, col in table_spec:
        if col is None or col not in bl_data.columns:
            print(f"{var_label:<32}{'--':<22}{'--':<22}{'--':<10}"); continue
        e1 = bl_data[bl_data[group_col]==1][col].dropna()
        e0 = bl_data[bl_data[group_col]==0][col].dropna()
        if len(e1)<3 or len(e0)<3:
            print(f"{var_label:<32}{'insufficient n':<44}"); continue
        _, p = stats.ttest_ind(e1, e0, equal_var=False)
        e1_str = f"{e1.mean():.1f} \u00b1 {e1.std():.1f}"
        e0_str = f"{e0.mean():.1f} \u00b1 {e0.std():.1f}"
        p_str = f"{p:.4f}{sig(p)}"
        print(f"{var_label:<32}{e1_str:<22}{e0_str:<22}{p_str:<10}")
        table_rows.append({'Variable':var_label,'Early_mean_SD':e1_str,'NonEarly_mean_SD':e0_str,'p':p,
                           'n_early':len(e1),'n_late':len(e0)})

    if 'Sex' in bl_data.columns:
        try:
            _sex_ct = pd.crosstab(bl_data[group_col], bl_data['Sex'])
            _, _p_sex, _, _ = stats.chi2_contingency(_sex_ct)
            e1_sex = bl_data[bl_data[group_col]==1]['Sex'].value_counts()
            e0_sex = bl_data[bl_data[group_col]==0]['Sex'].value_counts()
            e1_str = ':'.join(str(v) for v in e1_sex.values); e0_str = ':'.join(str(v) for v in e0_sex.values)
            p_str = f"{_p_sex:.4f}{sig(_p_sex)}"
            print(f"{'Sex (M:F)':<32}{e1_str:<22}{e0_str:<22}{p_str:<10}")
            table_rows.append({'Variable':'Sex (M:F)','Early_mean_SD':e1_str,'NonEarly_mean_SD':e0_str,'p':_p_sex,
                               'n_early':int((bl_data[group_col]==1).sum()),'n_late':int((bl_data[group_col]==0).sum())})
        except Exception as e:
            print(f"  Sex chi-square failed: {e}")

    out_df = pd.DataFrame(table_rows)
    out_df.to_csv(f'{OUT_DIR}/{label}_9H_demographics.csv', index=False)
    out_df.to_excel(f'{OUT_DIR}/{label}_9H_demographics_formatted.xlsx', index=False)
    return out_df

vh_demo_9h = _run_demographics_9h(vh_bl_data_9d, 'early_VH', 'VH')
ci_demo_9h = _run_demographics_9h(ci_bl_data_9d, 'early_CI', 'CI')
