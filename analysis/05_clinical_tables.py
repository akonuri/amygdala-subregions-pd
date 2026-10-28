"""Demographic and clinical characteristics (Table 1)."""

DEMO_VARS = [('Age','Age'),('EDUCYRS','Education (years)'),('disease_duration_at_MRI_yrs','Disease duration (years)'),
             ('MoCA_total','MoCA'),('TrailsB_sec','TMT-B (sec)'),('UPDRS3_total','MDS-UPDRS III'),
             ('RBDSQ_total','RBDSQ'),('SCOPA_total','SCOPA-Autonomic'),('GDS_total','GDS'),
             ('ESS_total','ESS'),('UPSIT_total','UPSIT'),('STAI_state','STAI-State')]

demo_table_rows = []
for col, label in [(None,'N'),('Sex','Sex (M:F)')] + DEMO_VARS:
    row = {'Variable':label}
    for grp in GROUPS + ['All_PD']:
        sub = get_grp(bl_all, grp)
        if col is None: v = str(len(sub))
        elif col=='Sex': v = f"{(sub['Sex']=='M').sum()}:{(sub['Sex']=='F').sum()}"
        else:
            vv = sub[col].dropna() if col in sub.columns else pd.Series(dtype=float)
            v = f'{vv.mean():.1f}+/-{vv.std():.1f}' if len(vv)>0 else 'NA'
        row[grp] = v
    demo_table_rows.append(row)
demo_df = pd.DataFrame(demo_table_rows)
demo_df.to_csv(f'{OUT_DIR}/demographics_table.csv', index=False)
print(demo_df.to_string(index=False))

demo_sig_rows = []; kw_ps = []
for col, label in DEMO_VARS:
    grps = [get_grp(bl_all, st)[col].dropna().values for st in SUBTYPES]
    valid = [g for g in grps if len(g)>3]
    kw_ps.append(stats.kruskal(*valid)[1] if len(valid)>=2 else np.nan)
valid_mask = [not np.isnan(p) for p in kw_ps]
kw_p_fdr = [np.nan]*len(kw_ps)
if sum(valid_mask) > 0:
    vp = [p for p,m in zip(kw_ps, valid_mask) if m]
    _, fc, _, _ = multipletests(vp, method='fdr_bh')
    j = 0
    for i,m in enumerate(valid_mask):
        if m: kw_p_fdr[i] = fc[j]; j += 1
for (col, label), p_kw, p_fdr in zip(DEMO_VARS, kw_ps, kw_p_fdr):
    if np.isnan(p_kw): continue
    print(f"\n{label}: KW p={p_kw:.4f} p_fdr={p_fdr:.4f}{sig(p_fdr)}")
    grps_dict = {st: get_grp(bl_all, st)[col].dropna().values for st in SUBTYPES}
    all_v = np.concatenate(list(grps_dict.values()))
    all_l = sum([[st]*len(v) for st,v in grps_dict.items()], [])
    try:
        ph = posthoc_conover(pd.DataFrame({'val':all_v,'grp':all_l}), val_col='val', group_col='grp', p_adjust='fdr_bh')
        for s1,s2 in [('I_Mild','II_Intermediate'),('I_Mild','III_Diffuse'),('II_Intermediate','III_Diffuse')]:
            if s1 in ph.index and s2 in ph.columns:
                p2 = ph.loc[s1,s2]
                print(f"  {PL[s1]} vs {PL[s2]}: p_fdr={p2:.4f}{sig(p2)}")
                demo_sig_rows.append({'variable':label,'KW_p':p_kw,'KW_p_fdr':p_fdr,'comparison':f'{PL[s1]}_vs_{PL[s2]}','conover_p_fdr':p2})
    except: pass
pd.DataFrame(demo_sig_rows).to_csv(f'{OUT_DIR}/demographics_within_subtype_significance.csv', index=False)

PDHC_VARS = [('Age','Age'),('MoCA_total','MoCA'),('TrailsB_sec','TMT-B (sec)'),
             ('UPDRS3_total','MDS-UPDRS III'),('UPDRS2_total','MDS-UPDRS II'),
             ('RBDSQ_total','RBDSQ'),('SCOPA_total','SCOPA-Autonomic'),
             ('GDS_total','GDS'),('ESS_total','ESS'),('UPSIT_total','UPSIT'),('STAI_state','STAI-State')]

pdhc_rows = []; pdhc_ps = []; pdhc_meta = []
for col, label in PDHC_VARS:
    hc_v = get_grp(bl_all,'HC')[col].dropna().values
    pd_v = get_grp(bl_all,'All_PD')[col].dropna().values
    if len(hc_v)<3 or len(pd_v)<3: continue
    stat, p = stats.mannwhitneyu(hc_v, pd_v, alternative='two-sided')
    pdhc_ps.append(p); pdhc_meta.append((label, hc_v, pd_v, stat))
pdhc_p_fdr = multipletests(pdhc_ps, method='fdr_bh')[1] if pdhc_ps else []
for (label, hc_v, pd_v, stat), p, p_fdr in zip(pdhc_meta, pdhc_ps, pdhc_p_fdr):
    print(f"  {label:<16}: HC={hc_v.mean():.2f} | PD={pd_v.mean():.2f} | p={p:.4f} p_fdr={p_fdr:.4f}{sig(p_fdr)}")
    pdhc_rows.append({'variable':label,'HC_mean':hc_v.mean(),'PD_mean':pd_v.mean(),'p':p,'p_fdr':p_fdr,'n_HC':len(hc_v),'n_PD':len(pd_v)})
hc_sex = get_grp(bl_all,'HC')['Sex'].dropna(); pd_sex = get_grp(bl_all,'All_PD')['Sex'].dropna()
if len(hc_sex)>0 and len(pd_sex)>0:
    ct_sex = pd.DataFrame({'M':[hc_sex.eq('M').sum(),pd_sex.eq('M').sum()],'F':[hc_sex.eq('F').sum(),pd_sex.eq('F').sum()]}, index=['HC','PD'])
    _, p_sex, _, _ = stats.chi2_contingency(ct_sex)
    print(f"  {'Sex':<16}: p={p_sex:.4f}{sig(p_sex)}")
    pdhc_rows.append({'variable':'Sex','p':p_sex,'p_fdr':p_sex})
pd.DataFrame(pdhc_rows).to_csv(f'{OUT_DIR}/PD_vs_HC_clinical.csv', index=False)
