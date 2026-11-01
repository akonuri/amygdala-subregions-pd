"""Exploratory analysis: incident depression (GDS-15 >=5) and anxiety (STAI state+trait >=54), Cox models with delayed entry, unadjusted and adjusted for age and sex."""

def _find_col_15(df_, cands):
    for c in cands:
        if c in df_.columns: return c
    return None

_GDS_C, _STAI_T, _STAI_S, _STAI_TR = ['GDS_total', 'GDSTOT', 'GDS_TOTAL'], ['STAI_total', 'STAITOT', 'STAI_TOTAL'], ['STAI_state', 'STAI_STATE'], ['STAI_trait', 'STAI_TRAIT']
_base_15 = _merged_dx_fu[['PATNO', 'EVENT_ID', 'Visit', 'years_from_diagnosis']].copy()
_base_15['PATNO'] = _base_15['PATNO'].astype(str)

def _source_15():
    for name, src, keys in [('PPMI_mega_clinical', mega, ['PATNO', 'EVENT_ID']), ('PPMI_volumes_clinical', vol, ['PATNO', 'Visit'])]:
        g = _find_col_15(src, _GDS_C); st = _find_col_15(src, _STAI_T)
        ss = _find_col_15(src, _STAI_S); sr = _find_col_15(src, _STAI_TR)
        if g or st or (ss and sr):
            cols = [c for c in [g, st, ss, sr] if c]
            d = src[keys + cols].copy(); d['PATNO'] = d['PATNO'].astype(str)
            for c in cols: d[c] = pd.to_numeric(d[c], errors='coerce')
            d = d.groupby(keys, as_index=False).first()
            print(f"Mood/anxiety source: {name} | GDS={g} STAI_total={st} STAI_state={ss} STAI_trait={sr}")
            return _base_15.merge(d, on=keys, how='left'), g, st, ss, sr
    return None, None, None, None, None

_mood_dx_15, _c_gds, _c_stai, _c_stais, _c_stait = _source_15()

def _spec_15(col, thr, total=None):
    def f():
        d = _mood_dx_15[['PATNO', 'years_from_diagnosis']].copy()
        v = (_mood_dx_15[total[0]] + _mood_dx_15[total[1]]) if total else _mood_dx_15[col]
        d['pos'] = np.where(v.notna(), (v >= thr).astype(float), np.nan)
        return d
    return f

_defs_15 = []
if _mood_dx_15 is not None:
    if _c_gds: _defs_15.append(('DEP_GDS5', 'Depression, GDS-15 >=5', _spec_15(_c_gds, 5)))
    if _c_stai: _defs_15.append(('ANX_STAI54', 'Anxiety, STAI total >=54', _spec_15(_c_stai, 54)))
    elif _c_stais and _c_stait: _defs_15.append(('ANX_STAI54', 'Anxiety, STAI state+trait >=54', _spec_15(None, 54, total=(_c_stais, _c_stait))))

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
    print("GDS or STAI not found in the clinical files, or too few events; nothing run.")
