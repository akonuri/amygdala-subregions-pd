"""Amygdala volumes across clinical subtypes (Table 2)."""

print("\nCross-sectional subtype comparison I vs II vs III (sectors, FDR family of 3)")
sub_rows = []; sub_ps = []; sub_meta = []
for col, label in SECTORS:
    grps = [get_grp(bl_all, st)[col].dropna().values for st in SUBTYPES]
    valid = [g for g in grps if len(g)>3]
    if len(valid)<2: continue
    _, p_kw = stats.kruskal(*valid)
    sub_ps.append(p_kw); sub_meta.append((col,label,grps))
sub_p_fdr = multipletests(sub_ps, method='fdr_bh')[1] if sub_ps else []
sub_pairwise = []
for (col,label,grps), p_kw, p_fdr in zip(sub_meta, sub_ps, sub_p_fdr):
    means = [f'{g.mean():.1f}' if len(g)>0 else 'NA' for g in grps]
    print(f"  {label:<8}: I={means[0]} | II={means[1]} | III={means[2]} mm3  KW p={p_kw:.4f} p_fdr={p_fdr:.4f}{sig(p_fdr)}")
    all_v = np.concatenate(grps); all_l = sum([[st]*len(v) for st,v in zip(SUBTYPES,grps)],[])
    try:
        ph = posthoc_conover(pd.DataFrame({'val':all_v,'grp':all_l}), val_col='val', group_col='grp', p_adjust='fdr_bh')
        for s1,s2 in [('I_Mild','II_Intermediate'),('I_Mild','III_Diffuse'),('II_Intermediate','III_Diffuse')]:
            if s1 in ph.index and s2 in ph.columns:
                p2 = ph.loc[s1,s2]
                if p2 < 0.10: print(f"    {PL[s1]} vs {PL[s2]}: p_fdr={p2:.4f}{sig(p2)}")
                sub_pairwise.append({'region':label,'col':col,'comparison':f'{PL[s1]}_vs_{PL[s2]}','conover_p_fdr':p2})
    except: pass
    sub_rows.append({'region':label,'col':col,'KW_p':p_kw,'KW_p_fdr':p_fdr})
pd.DataFrame(sub_rows).to_csv(f'{OUT_DIR}/VH11b_subtype_sectors.csv', index=False)
pd.DataFrame(sub_pairwise).to_csv(f'{OUT_DIR}/VH11b_subtype_sectors_pairwise.csv', index=False)

any_subtype_sec_sig = (pd.DataFrame(sub_rows)['KW_p_fdr'] < 0.05).any() if sub_rows else False
if any_subtype_sec_sig:
    print("\n  Post-hoc: nucleus-level subtype comparison (own FDR family of 9)")
    nuc_sub_rows = []; nuc_sub_ps = []; nuc_sub_meta = []
    for col, label in NUCLEI:
        grps = [get_grp(bl_all, st)[col].dropna().values for st in SUBTYPES]
        valid = [g for g in grps if len(g)>3]
        if len(valid)<2: continue
        _, p_kw = stats.kruskal(*valid)
        nuc_sub_ps.append(p_kw); nuc_sub_meta.append((col,label,grps))
    nuc_sub_p_fdr = multipletests(nuc_sub_ps, method='fdr_bh')[1] if nuc_sub_ps else []
    for (col,label,grps), p_kw, p_fdr in zip(nuc_sub_meta, nuc_sub_ps, nuc_sub_p_fdr):
        means = [f'{g.mean():.1f}' if len(g)>0 else 'NA' for g in grps]
        print(f"    {label:<28}: I={means[0]} | II={means[1]} | III={means[2]} mm3  KW p={p_kw:.4f} p_fdr={p_fdr:.4f}{sig(p_fdr)}")
        nuc_sub_rows.append({'region':label,'col':col,'KW_p':p_kw,'KW_p_fdr':p_fdr})
    pd.DataFrame(nuc_sub_rows).to_csv(f'{OUT_DIR}/VH11b_subtype_nuclei_posthoc.csv', index=False)
else:
    print("\n  No sector cleared FDR significance for subtype -- nucleus-level post-hoc skipped.")

print("\nLongitudinal sector x subtype x time interaction (sectors, FDR family of 3)")
def lme_sector_subtype_time(data, col, time_col='years_from_first_mri'):
    sub = data[(data['Group']=='PD') & (data['Visit'].isin(VISITS))].copy()
    sub = sub.dropna(subset=[col, time_col, 'subtype']).copy()
    sub = sub[sub['subtype'].isin(SUBTYPES)].copy()
    multi = sub.groupby('PATNO').filter(lambda x: len(x) > 1)
    if multi['PATNO'].nunique() < 15: return None
    col_s = col.replace('-','_')
    multi = multi.copy(); multi[col_s] = multi[col]
    multi['subtype_cat'] = pd.Categorical(multi['subtype'], categories=SUBTYPES)
    try:
        mdl = smf.mixedlm(f'{col_s} ~ {time_col} * C(subtype_cat, Treatment("I_Mild"))',
                           data=multi, groups=multi['PATNO']).fit(reml=True)
        p = mdl.params; pv = mdl.pvalues
        ii_key = [k for k in p.index if 'II_Intermediate' in k and time_col in k]
        iii_key = [k for k in p.index if 'III_Diffuse' in k and time_col in k]
        return {'n':multi['PATNO'].nunique(),
                'slope_I':p.get(time_col,np.nan),'p_slope_I':pv.get(time_col,np.nan),
                'interaction_II':p[ii_key[0]] if ii_key else np.nan,'p_interaction_II':pv[ii_key[0]] if ii_key else np.nan,
                'interaction_III':p[iii_key[0]] if iii_key else np.nan,'p_interaction_III':pv[iii_key[0]] if iii_key else np.nan}
    except: return None

def _omnibus_subtype_time_11c(data, col, time_col='years_from_first_mri'):
    sub = data[(data['Group']=='PD') & (data['Visit'].isin(VISITS))].copy()
    sub = sub.dropna(subset=[col, time_col, 'subtype']).copy()
    sub = sub[sub['subtype'].isin(SUBTYPES)].copy()
    multi = sub.groupby('PATNO').filter(lambda x: len(x) > 1).copy()
    if multi['PATNO'].nunique() < 15: return np.nan, np.nan
    col_s = col.replace('-','_'); multi[col_s] = multi[col]
    multi['subtype_cat'] = pd.Categorical(multi['subtype'], categories=SUBTYPES)
    try:
        full = smf.mixedlm(f'{col_s} ~ {time_col} * C(subtype_cat, Treatment("I_Mild"))', data=multi, groups=multi['PATNO']).fit(reml=False)
        red  = smf.mixedlm(f'{col_s} ~ {time_col} + C(subtype_cat, Treatment("I_Mild"))', data=multi, groups=multi['PATNO']).fit(reml=False)
        lr = max(2 * (full.llf - red.llf), 0.0)
        return lr, stats.chi2.sf(lr, len(SUBTYPES) - 1)
    except Exception:
        return np.nan, np.nan

long_sub_rows = []
for col, label in SECTORS:
    res = lme_sector_subtype_time(df, col)
    if res is None: continue
    lr_, po_ = _omnibus_subtype_time_11c(df, col)
    long_sub_rows.append({'sector':label,'col':col, 'omnibus_chi2':lr_, 'p_omnibus':po_, **res})
_ls = pd.DataFrame(long_sub_rows)
if len(_ls):
    _ls['p_omnibus_fdr'] = multipletests(_ls['p_omnibus'].fillna(1.0).values, method='fdr_bh')[1]
    _ls['p_fdr_I_vs_II'] = np.nan; _ls['p_fdr_I_vs_III'] = np.nan
    for _i in _ls.index:
        if _ls.loc[_i, 'p_omnibus_fdr'] < 0.05:
            _pf = multipletests([_ls.loc[_i, 'p_interaction_II'], _ls.loc[_i, 'p_interaction_III']], method='fdr_bh')[1]
            _ls.loc[_i, 'p_fdr_I_vs_II'] = _pf[0]; _ls.loc[_i, 'p_fdr_I_vs_III'] = _pf[1]
    for _r in _ls.itertuples():
        print(f"  {_r.sector:<8} (n={_r.n}): omnibus chi2({len(SUBTYPES)-1})={_r.omnibus_chi2:.2f} p={_r.p_omnibus:.4f} p_fdr={_r.p_omnibus_fdr:.4f}{sig(_r.p_omnibus_fdr)} | I(ref) slope={_r.slope_I:+.3f}")
        if _r.p_omnibus_fdr < 0.05:
            print(f"      I_vs_II={_r.interaction_II:+.3f} p_fdr={_r.p_fdr_I_vs_II:.4f}{sig(_r.p_fdr_I_vs_II)} | I_vs_III={_r.interaction_III:+.3f} p_fdr={_r.p_fdr_I_vs_III:.4f}{sig(_r.p_fdr_I_vs_III)} (family of 2)")
        else:
            print("      pairwise contrasts not interpreted (omnibus not FDR-significant)")
_ls.to_csv(f'{OUT_DIR}/VH11c_longitudinal_subtype_time.csv', index=False)
