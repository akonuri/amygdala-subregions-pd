"""Group-level source data for Figures 1 to 3 (for the Zenodo data record). Writes no participant-level values."""

_REGIONS_13 = [(whole_col, whole_label)] + SECTORS

def _raincloud_summary_13(bl_data, group_col, label):
    rows = []
    for col, wl in _REGIONS_13:
        for g, gname in [(0, f'Non-Early {label}'), (1, f'Early {label}')]:
            v = bl_data.loc[bl_data[group_col] == g, col].dropna()
            if len(v) == 0: continue
            rows.append({'symptom': label, 'region': wl, 'group': gname, 'n': len(v),
                         'mean_mm3': v.mean(), 'sd_mm3': v.std(), 'median_mm3': v.median(),
                         'q1_mm3': v.quantile(0.25), 'q3_mm3': v.quantile(0.75)})
    return pd.DataFrame(rows)

def _emm_table_13(cache, label):
    rows = []
    for col, wl in _REGIONS_13:
        if col not in cache: continue
        c = cache[col]
        for i, yr in enumerate(c['years']):
            rows.append({'symptom': label, 'region': wl, 'year': yr,
                         'emm_nonearly_mm3': c['e0s'][i], 'se_nonearly': c['se0s'][i],
                         'emm_early_mm3': c['e1s'][i], 'se_early': c['se1s'][i],
                         'diff_early_minus_nonearly': c['e1s'][i] - c['e0s'][i],
                         'p_uncorrected': c['p_uncorrected'][i], 'p_fdr': c['p_fdr'][i]})
    return pd.DataFrame(rows)

def _km_table_13(label, grid=np.arange(0, 15.5, 0.5)):
    coh, _ = _incident_cohort_9g(_SURV_SPEC_9G[label](), label)
    rows = []
    for col, wl in SECTORS:
        s = coh.merge(bl_pd[['PATNO', col]].assign(PATNO=lambda x: x['PATNO'].astype(str)), on='PATNO').dropna(subset=[col])
        low = (s[col] <= s[col].median()).values
        for m, gname in [(low, 'Lower volume'), (~low, 'Higher volume')]:
            km = _KMF9().fit(s['duration'][m], s['event'][m], entry=s['entry'][m], label=gname)
            sf = km.survival_function_at_times(grid).values
            ci = km.confidence_interval_survival_function_.reindex(
                km.confidence_interval_survival_function_.index.union(grid)).ffill().loc[grid]
            ent = s['entry'][m].values; dur = s['duration'][m].values
            for t, sv, lo, hi in zip(grid, sf, ci.iloc[:, 0].values, ci.iloc[:, 1].values):
                rows.append({'symptom': label, 'region': wl, 'group': gname,
                             'years_from_diagnosis': t, 'survival': sv, 'ci_lower': lo, 'ci_upper': hi,
                             'n_at_risk': int(((ent <= t) & (dur >= t)).sum())})
    return pd.DataFrame(rows)

_src13 = {
    'fig1A_VH_baseline_summary.csv': _raincloud_summary_13(vh_bl_data_9d, 'early_VH', 'VH'),
    'fig1B_VH_EMM.csv': _emm_table_13(vh_emm_cache_9e, 'VH'),
    'fig2A_CI_baseline_summary.csv': _raincloud_summary_13(ci_bl_data_9d, 'early_CI', 'CI'),
    'fig2B_CI_EMM.csv': _emm_table_13(ci_emm_cache_9e, 'CI'),
    'fig3A_VH_KaplanMeier.csv': _km_table_13('VH'),
    'fig3B_CI_KaplanMeier.csv': _km_table_13('CI'),
}
_src_dir13 = os.path.join(OUT_DIR, 'figure_source_data')
os.makedirs(_src_dir13, exist_ok=True)
for _fn, _df in _src13.items():
    _df.to_csv(os.path.join(_src_dir13, _fn), index=False)
    print(f"  Saved: figure_source_data/{_fn} ({len(_df)} rows)")

_xl13 = f'{OUT_DIR}/AMYGDALA_PD_MANUSCRIPT.xlsx'
if os.path.exists(_xl13):
    with pd.ExcelWriter(_xl13, engine='openpyxl', mode='a', if_sheet_exists='replace') as _w:
        for _fn, _df in _src13.items():
            _df.to_excel(_w, sheet_name=_fn.replace('.csv', '')[:31], index=False)
    print(f"  Added figure source data sheets to {_xl13}")
