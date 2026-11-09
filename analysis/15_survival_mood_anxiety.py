"""Exploratory analysis: incident depression (GDS-15 >=5) and anxiety (STAI-State >=40), scored from the PPMI item-level files at every clinical visit. Cox models with delayed entry, unadjusted and adjusted for age and sex."""
import glob as _glob15

def _find_file_15(pattern):
    hits = sorted(_glob15.glob(f'{BASE_DIR}/PPMI_Clinical_13Apr2026/**/{pattern}', recursive=True))
    return hits[0] if hits else None

def _read_15(path):
    d = safe_read_csv(path, low_memory=False)
    d.columns = [c.replace('-', '_') for c in d.columns]
    d['PATNO'] = d['PATNO'].astype(str)
    return d

_GDS_POS = ['GDSSATIS', 'GDSGSPIR', 'GDSHAPPY', 'GDSALIVE', 'GDSENRGY']
_GDS_NEG = ['GDSDROPD', 'GDSEMPTY', 'GDSBORED', 'GDSAFRAD', 'GDSHLPLS', 'GDSHOME', 'GDSMEMRY', 'GDSWRTLS', 'GDSHOPLS', 'GDSBETER']

def _score_gds_15(d):
    items = [c for c in _GDS_POS + _GDS_NEG if c in d.columns]
    if len(items) < 15: return None
    x = d[items].apply(pd.to_numeric, errors='coerce').where(lambda v: v.isin([0, 1]))
    pts = pd.concat([1 - x[_GDS_POS], x[_GDS_NEG]], axis=1)
    n_ans = pts.notna().sum(axis=1)
    return (pts.sum(axis=1) * 15 / n_ans).where(n_ans >= 12)

_STAI_REV = [1, 2, 5, 8, 10, 11, 15, 16, 19, 20]

def _score_stai_state_15(d):
    items = [f'STAIAD{i}' for i in range(1, 21)]
    if not all(c in d.columns for c in items): return None
    x = d[items].apply(pd.to_numeric, errors='coerce').where(lambda v: v.isin([1, 2, 3, 4]))
    for i in _STAI_REV: x[f'STAIAD{i}'] = 5 - x[f'STAIAD{i}']
    n_ans = x.notna().sum(axis=1)
    return (x.sum(axis=1) * 20 / n_ans).where(n_ans >= 18)

_scores_15 = {}
_f_gds = _find_file_15('Geriatric_Depression_Scale*.csv')
_f_stai = _find_file_15('State*Trait_Anxiety*.csv')
print(f"GDS file: {_f_gds}\nSTAI file: {_f_stai}")
if _f_gds:
    _g = _read_15(_f_gds); _gs = _score_gds_15(_g)
    if _gs is not None and _gs.notna().any():
        _g['GDS15'] = _gs
        _scores_15['GDS15'] = _g[['PATNO', 'EVENT_ID', 'GDS15']].dropna().drop_duplicates(['PATNO', 'EVENT_ID'])
        print(f"  GDS-15 scored: {len(_scores_15['GDS15'])} visits, {_scores_15['GDS15']['PATNO'].nunique()} participants")
    else: print(f"  GDS items not found. Columns: {list(_g.columns)[:40]}")
if _f_stai:
    _st = _read_15(_f_stai); _sc = _score_stai_state_15(_st)
    if _sc is not None and _sc.notna().any():
        _st['STAI_S'] = _sc
        _scores_15['STAI_S'] = _st[['PATNO', 'EVENT_ID', 'STAI_S']].dropna().drop_duplicates(['PATNO', 'EVENT_ID'])
        print(f"  STAI-State scored: {len(_scores_15['STAI_S'])} visits, {_scores_15['STAI_S']['PATNO'].nunique()} participants")
    else: print(f"  STAI items not found. Columns: {list(_st.columns)[:40]}")

_base_15 = _merged_dx_fu[['PATNO', 'EVENT_ID', 'years_from_diagnosis']].copy()
_base_15['PATNO'] = _base_15['PATNO'].astype(str)
_mood_dx_15 = _base_15
for _k, _v in _scores_15.items():
    _mood_dx_15 = _mood_dx_15.merge(_v, on=['PATNO', 'EVENT_ID'], how='left')

def _spec_15(col, thr):
    def f():
        d = _mood_dx_15[['PATNO', 'years_from_diagnosis']].copy()
        v = _mood_dx_15[col]
        d['pos'] = np.where(v.notna(), (v >= thr).astype(float), np.nan)
        return d
    return f

_defs_15 = []
if 'GDS15' in _scores_15: _defs_15.append(('DEP_GDS5', 'Depression, GDS-15 >=5', _spec_15('GDS15', 5)))
if 'STAI_S' in _scores_15: _defs_15.append(('ANX_STAIS40', 'Anxiety, STAI-State >=40', _spec_15('STAI_S', 40)))

_res_15 = []
for _lab, _desc, _fn in _defs_15:
    print(f"\n{_lab} ({_desc})")
    _SURV_SPEC_9G[_lab] = _fn
    _chk, _ = _incident_cohort_9g(_fn(), _lab)
    if len(_chk) == 0 or _chk['event'].sum() < 10:
        print(f"  {_lab}: too few symptom-free participants or events, skipped"); continue
    try:
        _u = _run_survival_all_sectors_9g(None, _lab)
    except Exception as _e:
        print(f"  {_lab}: unadjusted survival failed: {_e!r}"); continue
    if not len(_u): continue
    _u = _u.assign(outcome=_lab, definition=_desc, model='unadjusted')
    _s = _chk.merge(bl_pd[['PATNO', 'Age', 'Sex'] + [c for c, _ in SECTORS]].assign(PATNO=lambda x: x['PATNO'].astype(str)),
                    on='PATNO', how='inner')
    _rows = []
    for _col, _wl in SECTORS:
        if _s.dropna(subset=[_col])['event'].sum() < 10: continue
        try:
            _r = _cox_adjusted_14(_s, _col); _r.update({'region': _wl}); _rows.append(_r)
        except Exception as _e:
            print(f"  {_wl}: adjusted model failed: {_e!r}")
    if _rows:
        _a = pd.DataFrame(_rows)
        _a['p_fdr'] = multipletests(_a['p'].values, method='fdr_bh')[1]
        print("  Adjusted for age and sex:")
        for r in _a.itertuples():
            print(f"    {r.region:<8}: HR/SD lower={r.hr_per_SD_lower:.3f} [{r.ci_lo:.3f}, {r.ci_hi:.3f}] p={r.p:.4f} p_fdr={r.p_fdr:.4f}{sig(r.p_fdr)}")
        _res_15.append(_a.assign(outcome=_lab, definition=_desc, model='age_sex_adjusted'))
    _res_15.append(_u)

if _res_15:
    _keep = ['outcome', 'definition', 'model', 'region', 'hr_per_SD_lower', 'ci_lo', 'ci_hi', 'p', 'p_fdr', 'n', 'events']
    _out15 = pd.concat([r[[c for c in _keep if c in r.columns]] for r in _res_15], ignore_index=True)
    _out15.to_csv(f'{OUT_DIR}/survival_mood_anxiety.csv', index=False)
    _xl15 = f'{OUT_DIR}/AMYGDALA_PD_MANUSCRIPT.xlsx'
    if os.path.exists(_xl15):
        with pd.ExcelWriter(_xl15, engine='openpyxl', mode='a', if_sheet_exists='replace') as _w:
            _out15.to_excel(_w, sheet_name='Survival_Mood_Anxiety', index=False)
    print("\nSaved: survival_mood_anxiety.csv")
else:
    print("No scorable GDS or STAI data found, or too few events; nothing run.")
