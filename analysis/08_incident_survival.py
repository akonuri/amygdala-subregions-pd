"""Incident VH and CI, Cox models with delayed entry and Kaplan-Meier analyses, plus VH threshold sensitivity analysis."""

from lifelines import KaplanMeierFitter as _KMF9, CoxPHFitter as _Cox9
from scipy import stats as _st9
from statsmodels.stats.multitest import multipletests as _mt9

def _logrank_delayed_entry_9g(entry, duration, event, group):
    entry = np.asarray(entry, float); duration = np.asarray(duration, float)
    event = np.asarray(event, int); group = np.asarray(group, int)
    O = E = V = 0.0
    for t in np.unique(duration[event == 1]):
        atrisk = (entry < t) & (duration >= t)
        n = atrisk.sum(); n1 = (atrisk & (group == 1)).sum()
        m = atrisk & (duration == t) & (event == 1)
        d = m.sum(); d1 = (m & (group == 1)).sum()
        if n < 2 or d == 0: continue
        O += d1; E += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    if V <= 0: return np.nan, np.nan
    chi = (O - E) ** 2 / V
    return chi, _st9.chi2.sf(chi, 1)

def _long_pos_vh_9g():
    d = _merged_dx_fu[['PATNO', 'years_from_diagnosis', 'U1_NP1HALL']].copy()
    v = pd.to_numeric(d['U1_NP1HALL'], errors='coerce')
    v = v.where(v != 101)
    d['pos'] = np.where(v.notna(), (v >= 1).astype(float), np.nan)
    return d[['PATNO', 'years_from_diagnosis', 'pos']]

_CI_RULES_9G = [('HVLT_DVT_TOTAL_RECALL', 35), ('HVLT_DVT_RECOG_DISC_INDEX', 35),
                ('SEMFLU_DVS_SFTANIM', 5.5), ('SDMT_DVT_SDM', 35),
                ('LNS_DVS_LNS', 5.5), ('DVS_JLO_MSSA', 5.5)]

def _long_pos_ci_9g():
    global _merged_dx_fu
    if 'DVS_JLO_MSSA' not in _merged_dx_fu.columns:
        _j = safe_read_csv(f'{BASE_DIR}/PPMI_Clinical_13Apr2026/Non-motor_Assessments/Benton_Judgement_of_Line_Orientation_12Apr2026.csv', low_memory=False)
        _j.columns = [c.replace('-', '_') for c in _j.columns]
        _j['PATNO'] = _j['PATNO'].astype(str)
        _j['DVS_JLO_MSSA'] = pd.to_numeric(_j['DVS_JLO_MSSA'], errors='coerce')
        _j = _j.dropna(subset=['DVS_JLO_MSSA']).drop_duplicates(['PATNO', 'EVENT_ID'])[['PATNO', 'EVENT_ID', 'DVS_JLO_MSSA']]
        _tmp = _merged_dx_fu.copy(); _tmp['PATNO'] = _tmp['PATNO'].astype(str)
        _merged_dx_fu = _tmp.merge(_j, on=['PATNO', 'EVENT_ID'], how='left')
        print(f"  JLO merged onto _merged_dx_fu: {_merged_dx_fu['DVS_JLO_MSSA'].notna().sum()} visits with a score")
    have = [(c, th) for c, th in _CI_RULES_9G if c in _merged_dx_fu.columns]
    print("  CI tests found in _merged_dx_fu:", [c for c, _ in have])
    if len(have) < 4:
        raise KeyError("fewer than 4 CI test columns on _merged_dx_fu")
    d = _merged_dx_fu[['PATNO', 'years_from_diagnosis'] + [c for c, _ in have]].copy()
    imp = pd.DataFrame(index=d.index)
    for c, th in have:
        x = pd.to_numeric(d[c], errors='coerce')
        imp[c] = np.where(x.notna(), (x <= th).astype(float), np.nan)
    n_av = imp.notna().sum(axis=1); n_im = imp.sum(axis=1)
    d['pos'] = np.where(n_av >= 4, (n_im >= 2).astype(float), np.nan)
    return d[['PATNO', 'years_from_diagnosis', 'pos']]

_SURV_SPEC_9G = {'VH': _long_pos_vh_9g, 'CI': _long_pos_ci_9g}

def _incident_cohort_9g(lp, label=''):
    lp = lp.copy()
    lp['PATNO'] = lp['PATNO'].astype(str)
    lp = lp[lp['PATNO'].isin(set(bl_pd['PATNO'].astype(str)))]
    lp = lp.dropna(subset=['years_from_diagnosis', 'pos'])
    lp = lp.drop_duplicates(['PATNO', 'years_from_diagnosis', 'pos'])
    lp = lp.sort_values(['PATNO', 'years_from_diagnosis'])
    rows = []; flow = {'with_valid_assessment': 0, 'prevalent_at_first': 0,
                       'single_assessment': 0, 'included': 0, 'events': 0}
    for pid, g in lp.groupby('PATNO'):
        t = g['years_from_diagnosis'].values; p = g['pos'].values
        flow['with_valid_assessment'] += 1
        if p[0] >= 1:
            flow['prevalent_at_first'] += 1; continue
        k = np.where(p >= 1)[0]
        if len(k): dur, ev = t[k[0]], 1
        else:      dur, ev = t[-1], 0
        if dur <= t[0]:
            flow['single_assessment'] += 1; continue
        rows.append({'PATNO': pid, 'entry': t[0], 'duration': dur, 'event': ev})
        flow['included'] += 1; flow['events'] += ev
    coh = pd.DataFrame(rows)
    flow['entry_gt_5y'] = int((coh['entry'] > 5).sum()) if len(coh) else 0
    flow['entry_gt_10y'] = int((coh['entry'] > 10).sum()) if len(coh) else 0
    print(f"  [{label}] cohort flow: {flow}")
    if len(coh):
        print(f"  [{label}] entry (yrs from dx) median={coh['entry'].median():.2f} "
              f"range=[{coh['entry'].min():.2f}, {coh['entry'].max():.2f}] | entry<0: {(coh['entry']<0).sum()}")
    return coh, flow

def _run_cox_entry_9g(surv_df, region_col):
    tmp = surv_df.dropna(subset=[region_col, 'duration', 'event']).copy()
    tmp['region_z'] = (tmp[region_col] - tmp[region_col].mean()) / tmp[region_col].std()
    cols = ['duration', 'event', 'region_z'] + (['entry'] if 'entry' in tmp.columns else [])
    cph = _Cox9()
    if 'entry' in tmp.columns:
        cph.fit(tmp[cols], duration_col='duration', event_col='event', entry_col='entry')
    else:
        cph.fit(tmp[cols], duration_col='duration', event_col='event')
    hr = np.exp(cph.params_['region_z'])
    lo, hi = np.exp(cph.confidence_intervals_.loc['region_z'])
    p = cph.summary.loc['region_z', 'p']
    p_ph = np.nan
    from lifelines.statistics import proportional_hazard_test
    try:
        p_ph = proportional_hazard_test(cph, tmp[cols], time_transform='rank').summary['p'].iloc[0]
    except Exception as e1:
        try:
            _c2cols = ['duration', 'event', 'region_z']
            _c2 = _Cox9().fit(tmp[_c2cols], duration_col='duration', event_col='event')
            p_ph = proportional_hazard_test(_c2, tmp[_c2cols], time_transform='rank').summary['p'].iloc[0]
        except Exception as e2:
            print(f"    (PH test failed, {type(e2).__name__}: {str(e2)[:90]})")
    return 1 / hr, 1 / hi, 1 / lo, p, len(tmp), int(tmp['event'].sum()), p_ph

def _run_survival_all_sectors_9g(groups_df, label):
    spec = _SURV_SPEC_9G.get(label)
    if spec is None:
        raise KeyError(f"no incident-cohort spec registered for {label}")
    print(f"{label}: incident survival, symptom-free at first assessment, left-truncated at first assessment")
    coh, flow = _incident_cohort_9g(spec(), label)
    pd.DataFrame([flow]).to_csv(f'{OUT_DIR}/{label}_9G_cohort_flow.csv', index=False)
    rows = []
    for col, wlabel in SECTORS:
        s = coh.merge(bl_pd[['PATNO', col]].assign(PATNO=lambda x: x['PATNO'].astype(str)), on='PATNO', how='inner')
        s = s.dropna(subset=[col])
        if s['event'].sum() < 10:
            print(f"  {wlabel}: fewer than 10 events, skipped"); continue
        hr, lo, hi, p, n, ne, p_ph = _run_cox_entry_9g(s, col)
        low = (s[col] <= s[col].median()).astype(int).values
        _, p_lr = _logrank_delayed_entry_9g(s['entry'], s['duration'], s['event'], low)
        med = {}
        for gname, gv in [('low', 1), ('high', 0)]:
            m = low == gv
            km = _KMF9().fit(s['duration'][m], s['event'][m], entry=s['entry'][m])
            med[gname] = km.median_survival_time_
        rows.append({'region': wlabel, 'col': col, 'hr_per_SD_lower': hr, 'ci_lo': lo, 'ci_hi': hi,
                     'p': p, 'logrank_p': p_lr, 'median_low': med['low'], 'median_high': med['high'],
                     'n': n, 'events': ne, 'p_ph': p_ph})
    res = pd.DataFrame(rows)
    if len(res):
        res['p_fdr'] = _mt9(res['p'].values, method='fdr_bh')[1]
        res['logrank_p_fdr'] = _mt9(res['logrank_p'].fillna(1.0).values, method='fdr_bh')[1]
        for r in res.itertuples():
            print(f"  {r.region:<8}: HR/SD lower={r.hr_per_SD_lower:.3f} [{r.ci_lo:.3f}, {r.ci_hi:.3f}] "
                  f"p={r.p:.4f} p_fdr={r.p_fdr:.4f}{sig(r.p_fdr)} | log-rank(delayed entry) p={r.logrank_p:.4f} p_fdr={r.logrank_p_fdr:.4f}{sig(r.logrank_p_fdr)} "
                  f"| median low={r.median_low:.1f} high={r.median_high:.1f} | n={r.n} events={r.events} | PH p={r.p_ph:.3f}")
    res.to_csv(f'{OUT_DIR}/{label}_9G_survival.csv', index=False)
    return res

vh_surv_results_9g = _run_survival_all_sectors_9g(vh_early_groups, 'VH')
ci_surv_results_9g = _run_survival_all_sectors_9g(ci_early_groups, 'CI')

_vh_sens_results = {}
for _thr in [2]:
    _lab = f'VHge{_thr}'
    print(f"\n{_lab}: NP1HALL >= {_thr}")
    _d = _merged_dx_fu.dropna(subset=['U1_NP1HALL']).copy()
    _d['U1_NP1HALL'] = pd.to_numeric(_d['U1_NP1HALL'], errors='coerce')
    _d = _d[_d['U1_NP1HALL'].notna() & (_d['U1_NP1HALL'] != 101)].copy()
    _d['outcome_pos'] = (_d['U1_NP1HALL'] >= _thr).astype(float)
    _grp = _build_early_groups(_d, 'outcome_pos', f'{_lab} (NP1HALL>={_thr})')
    _grp.to_excel(f'{OUT_DIR}/{_lab}_early_groups.xlsx', index=False)
    _gc = f'early_{_lab}'
    try:
        _bl, _sec, _nuc = _run_cross_sectional_9d(_grp, _gc, _lab)
        _run_longitudinal_9e(_grp, _gc, _lab)
        _run_moca_sensitivity_9f(_bl, _gc, _lab)
    except Exception as _e:
        print(f"  {_lab}: early-group analyses failed: {_e!r}")
    def _vh_spec(thr=_thr):
        d = _merged_dx_fu[['PATNO','years_from_diagnosis','U1_NP1HALL']].copy()
        v = pd.to_numeric(d['U1_NP1HALL'], errors='coerce')
        v = v.where(v != 101)
        d['pos'] = np.where(v.notna(), (v >= thr).astype(float), np.nan)
        return d[['PATNO','years_from_diagnosis','pos']]
    _SURV_SPEC_9G[_lab] = _vh_spec
    try:
        _vh_sens_results[_lab] = _run_survival_all_sectors_9g(_grp, _lab)
    except Exception as _e:
        print(f"  {_lab}: survival failed: {_e!r}")
