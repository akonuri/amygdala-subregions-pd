"""Sensitivity analysis: Cox models for incident VH and CI adjusted for age and sex at baseline MRI."""

def _cox_adjusted_14(s, col):
    t = s.dropna(subset=[col, 'Age', 'Sex', 'duration', 'event', 'entry']).copy()
    t['region_z'] = (t[col] - t[col].mean()) / t[col].std()
    t['age_10y'] = t['Age'] / 10.0
    t['male'] = (t['Sex'] == 'M').astype(float)
    cols = ['duration', 'event', 'entry', 'region_z', 'age_10y', 'male']
    cph = _Cox9().fit(t[cols], duration_col='duration', event_col='event', entry_col='entry')
    ci = np.exp(cph.confidence_intervals_)
    out = {'hr_per_SD_lower': 1 / np.exp(cph.params_['region_z']),
           'ci_lo': 1 / ci.loc['region_z'].iloc[1], 'ci_hi': 1 / ci.loc['region_z'].iloc[0],
           'p': cph.summary.loc['region_z', 'p'],
           'hr_age_per_10y': np.exp(cph.params_['age_10y']), 'p_age': cph.summary.loc['age_10y', 'p'],
           'hr_male': np.exp(cph.params_['male']), 'p_male': cph.summary.loc['male', 'p'],
           'n': len(t), 'events': int(t['event'].sum())}
    from lifelines.statistics import proportional_hazard_test
    try:
        out['p_ph_global'] = proportional_hazard_test(cph, t[cols], time_transform='rank').summary['p'].min()
    except Exception:
        try:
            _c = [c for c in cols if c != 'entry']
            _m = _Cox9().fit(t[_c], duration_col='duration', event_col='event')
            out['p_ph_global'] = proportional_hazard_test(_m, t[_c], time_transform='rank').summary['p'].min()
        except Exception:
            out['p_ph_global'] = np.nan
    return out

_adj_rows_14 = []
for _lab in [l for l in ['VH', 'CI', 'VHge2'] if l in _SURV_SPEC_9G]:
    print(f"\n{_lab}: Cox adjusted for age and sex (per SD lower volume)")
    _coh14, _ = _incident_cohort_9g(_SURV_SPEC_9G[_lab](), _lab)
    _bl14 = bl_pd[['PATNO', 'Age', 'Sex'] + [c for c, _ in SECTORS]].assign(PATNO=lambda x: x['PATNO'].astype(str))
    _s14 = _coh14.merge(_bl14, on='PATNO', how='inner')
    _rows = []
    for _col, _wl in SECTORS:
        if _s14.dropna(subset=[_col])['event'].sum() < 10: continue
        _r = _cox_adjusted_14(_s14, _col); _r.update({'outcome': _lab, 'region': _wl}); _rows.append(_r)
    if not _rows: continue
    _df = pd.DataFrame(_rows)
    _df['p_fdr'] = multipletests(_df['p'].values, method='fdr_bh')[1]
    for r in _df.itertuples():
        print(f"  {r.region:<8}: HR/SD lower={r.hr_per_SD_lower:.3f} [{r.ci_lo:.3f}, {r.ci_hi:.3f}] p={r.p:.4f} "
              f"p_fdr={r.p_fdr:.4f}{sig(r.p_fdr)} | age HR/10y={r.hr_age_per_10y:.2f} p={r.p_age:.4f} | "
              f"male HR={r.hr_male:.2f} p={r.p_male:.4f} | n={r.n} events={r.events} | PH p(min)={r.p_ph_global:.3f}")
    _adj_rows_14.append(_df)

if _adj_rows_14:
    _adj14 = pd.concat(_adj_rows_14, ignore_index=True)
    _adj14 = _adj14[['outcome', 'region', 'hr_per_SD_lower', 'ci_lo', 'ci_hi', 'p', 'p_fdr', 'hr_age_per_10y', 'p_age',
                     'hr_male', 'p_male', 'n', 'events', 'p_ph_global']]
    _adj14.to_csv(f'{OUT_DIR}/survival_age_sex_adjusted.csv', index=False)
    _xl14 = f'{OUT_DIR}/AMYGDALA_PD_MANUSCRIPT.xlsx'
    if os.path.exists(_xl14):
        with pd.ExcelWriter(_xl14, engine='openpyxl', mode='a', if_sheet_exists='replace') as _w:
            _adj14.to_excel(_w, sheet_name='Survival_AgeSex_Adjusted', index=False)
    print(f"\nSaved: survival_age_sex_adjusted.csv")
