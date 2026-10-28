"""Figures 1 to 3."""

_NOTEARLY_LINE = '#a8895f'; _NOTEARLY_FILL = '#e7d9cc'
_EARLY_LINE = '#3f7d73'; _EARLY_FILL = '#b9d9d5'

def _draw_raincloud_9eplot(ax, groups_data, colors, positions, width=0.32, point_size=10, point_alpha=0.28):
    from matplotlib.patches import Polygon
    for data, (fill, line), pos in zip(groups_data, colors, positions):
        data = np.asarray(data); data = data[~np.isnan(data)]
        if len(data) < 3: continue
        kde = stats.gaussian_kde(data)
        y_grid = np.linspace(data.min(), data.max(), 200)
        density = kde(y_grid); density = density/density.max()*width
        verts = list(zip(pos+density, y_grid)) + [(pos,y_grid[-1]), (pos,y_grid[0])]
        ax.add_patch(Polygon(verts, facecolor=fill, edgecolor=line, alpha=0.8, linewidth=1.2, zorder=2))
        jitter = np.random.RandomState(42).uniform(-width*0.55, -width*0.15, size=len(data))
        ax.scatter(pos+jitter, data, s=point_size, color=line, alpha=point_alpha, edgecolors='none', zorder=1)
        ax.boxplot([data], positions=[pos-width*0.62], widths=width*0.26, vert=True, patch_artist=True,
                  showfliers=False, zorder=3, boxprops=dict(facecolor=fill,edgecolor=line,linewidth=1.2),
                  medianprops=dict(color=line,linewidth=2.0), whiskerprops=dict(color=line,linewidth=1.1),
                  capprops=dict(color=line,linewidth=1.1))

def _add_sig_bracket_9eplot(ax, positions, marker, groups_data):
    if not marker: return
    all_vals = np.concatenate([np.asarray(g)[~np.isnan(np.asarray(g))] for g in groups_data if len(g)>0])
    if len(all_vals)==0: return
    y_top = all_vals.max(); y_range = all_vals.max()-all_vals.min()
    bracket_y = y_top+0.09*y_range; bracket_h = 0.025*y_range
    x_left, x_right = positions[0]-0.18, positions[-1]-0.18
    ax.plot([x_left,x_left,x_right,x_right],[bracket_y,bracket_y+bracket_h,bracket_y+bracket_h,bracket_y],
           color='#222222', linewidth=1.6, zorder=6, clip_on=False)
    ax.text((x_left+x_right)/2, bracket_y+bracket_h*1.4, marker, ha='center', va='bottom',
           fontsize=20, color='#222222', zorder=6, fontweight='bold')
    ax.set_ylim(top=bracket_y+bracket_h*6.5 if ax.get_ylim()[1] < bracket_y+bracket_h*6.5 else ax.get_ylim()[1])

def _build_raincloud_emm_figure(bl_data_9d, long_data_9e, emm_cache_9e, group_col, group_label, outpath, sec_results_9d=None, use_uncorrected_p=False):
    _sec_pfdr_lookup = {}
    if sec_results_9d:
        _sec_pfdr_lookup = {r['col']: r.get('p_fdr', np.nan) for r in sec_results_9d}
    _fig = plt.figure(figsize=(20, 11))
    _gsA = _fig.add_gridspec(2, 2, left=0.05, right=0.47, top=0.90, bottom=0.12, hspace=0.42, wspace=0.30)
    _gsB = _fig.add_gridspec(2, 2, left=0.55, right=0.98, top=0.90, bottom=0.12, hspace=0.42, wspace=0.30)
    _positions = np.array([0,1])
    _MAIN_REGIONS = [(whole_col, whole_label)] + SECTORS

    for _i, (_col, _label) in enumerate(_MAIN_REGIONS):
        _r, _c = divmod(_i, 2)
        _ax = _fig.add_subplot(_gsA[_r, _c]); _ax.set_box_aspect(1)
        _neg_d = bl_data_9d[bl_data_9d[group_col]==0][_col].dropna().values
        _pos_d = bl_data_9d[bl_data_9d[group_col]==1][_col].dropna().values
        if _col == whole_col:
            _, _p_cs = stats.ttest_ind(_pos_d, _neg_d, equal_var=False)
        else:
            _p_cs = _sec_pfdr_lookup.get(_col, np.nan)
        _sig_marker_cs = '***' if _p_cs<0.001 else '**' if _p_cs<0.01 else '*' if _p_cs<0.05 else ''
        _draw_raincloud_9eplot(_ax, [_neg_d, _pos_d], [(_NOTEARLY_FILL,_NOTEARLY_LINE),(_EARLY_FILL,_EARLY_LINE)], _positions)
        _add_sig_bracket_9eplot(_ax, _positions, _sig_marker_cs, [_neg_d, _pos_d])
        _ax.set_title(_label, fontsize=20, pad=18)
        _ax.set_xticks(_positions-0.18); _ax.set_xticklabels([f'Non-early {group_label}', f'Early {group_label}'], fontsize=15)
        _ax.set_ylabel('Volume (mm\u00b3)', fontsize=16)
        _ax.tick_params(axis='both', labelsize=13)
        _ax.set_xlim(-0.75, 1.55); _ax.grid(axis='y', alpha=0.12, linewidth=0.6)
        if _i == 0: _ax.text(-0.26, 1.14, 'A', transform=_ax.transAxes, fontsize=26, fontweight='bold')

    for _i, (_col, _label) in enumerate(_MAIN_REGIONS):
        _r, _c = divmod(_i, 2)
        _ax = _fig.add_subplot(_gsB[_r, _c]); _ax.set_box_aspect(1)
        if _col not in emm_cache_9e:
            _ax.set_title(f"{_label}\n(insufficient data)", fontsize=18); continue
        _cached = emm_cache_9e[_col]
        _e0s, _se0s, _e1s, _se1s = _cached['e0s'], _cached['se0s'], _cached['e1s'], _cached['se1s']
        _p_for_stars = _cached.get('p_uncorrected', _cached['p_fdr']) if use_uncorrected_p else _cached['p_fdr']

        _x_grid = np.array([0,1,2,3,4])
        _ax.plot(_x_grid, _e0s, '-', color=_NOTEARLY_LINE, linewidth=4.5, zorder=4, label=f'Non-early {group_label}')
        _ax.plot(_x_grid, _e0s, 'o', color=_NOTEARLY_LINE, markersize=11, markeredgecolor='white', markeredgewidth=1.4, zorder=5)
        _ax.plot(_x_grid, _e1s, '-', color=_EARLY_LINE, linewidth=4.5, zorder=4, label=f'Early {group_label}')
        _ax.plot(_x_grid, _e1s, 'o', color=_EARLY_LINE, markersize=11, markeredgecolor='white', markeredgewidth=1.4, zorder=5)

        for _t, _pf, _ee0, _ee1 in zip([0,1,2,3,4], _p_for_stars, _e0s, _e1s):
            if _pf < 0.05:
                _ax.annotate('*', xy=(_t, max(_ee0,_ee1)), xytext=(0,12), textcoords='offset points',
                           ha='center', fontsize=22, fontweight='bold', color='#222222', zorder=6)

        _all_vals_b = np.array(_e0s + _e1s)
        _y_range_b = _all_vals_b.max() - _all_vals_b.min()
        _ax.set_ylim(_all_vals_b.min() - 0.6*_y_range_b, _all_vals_b.max() + 0.6*_y_range_b)
        _ax.set_title(_label, fontsize=20)
        _ax.set_xlabel('Years from first MRI', fontsize=16)
        _ax.set_ylabel('Estimated marginal mean (mm\u00b3)', fontsize=16)
        _ax.tick_params(axis='both', labelsize=13)
        _ax.set_xticks([0,1,2,3,4]); _ax.grid(axis='y', alpha=0.12, linewidth=0.6)
        if _i == 0: _ax.text(-0.26, 1.14, 'B', transform=_ax.transAxes, fontsize=26, fontweight='bold')

    _handles, _labels = _ax.get_legend_handles_labels()
    _fig.legend(_handles, _labels, loc='lower center', ncol=2, fontsize=18, frameon=False,
               bbox_to_anchor=(0.5, 0.005), handlelength=2.5, columnspacing=3.0)
    plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {outpath}")

_build_raincloud_emm_figure(vh_bl_data_9d, vh_long_data_9e, vh_emm_cache_9e, 'early_VH', 'VH',
                            f'{PLOT_DIR}/fig_early_VH_raincloud_EMM_sidebyside.png', sec_results_9d=vh_sec_results_9d, use_uncorrected_p=False)
_build_raincloud_emm_figure(ci_bl_data_9d, ci_long_data_9e, ci_emm_cache_9e, 'early_CI', 'CI',
                            f'{PLOT_DIR}/fig_early_CI_raincloud_EMM_sidebyside.png', sec_results_9d=ci_sec_results_9d, use_uncorrected_p=False)

import matplotlib.pyplot as _plt9
_res9 = {'VH': vh_surv_results_9g, 'CI': ci_surv_results_9g}
_fig9, _axes9 = _plt9.subplots(2, 3, figsize=(15, 8.5), sharey=True)
for _i, (_lab, _full) in enumerate([('VH', 'VH'), ('CI', 'CI')]):
    _coh, _ = _incident_cohort_9g(_SURV_SPEC_9G[_lab](), _lab)
    for _j, (_col, _wl) in enumerate(SECTORS):
        _ax = _axes9[_i, _j]
        _s = _coh.merge(bl_pd[['PATNO', _col]].assign(PATNO=lambda x: x['PATNO'].astype(str)), on='PATNO').dropna(subset=[_col])
        _low = (_s[_col] <= _s[_col].median()).values
        for _m, _nm, _c in [(_low, 'Lower volume', '#D55E00'), (~_low, 'Higher volume', '#0072B2')]:
            _km = _KMF9().fit(_s['duration'][_m], _s['event'][_m], entry=_s['entry'][_m], label=_nm)
            _km.plot_survival_function(ax=_ax, ci_show=True, color=_c, linewidth=1.6)
        _r = _res9[_lab]
        _r = _r[_r['region'] == _wl] if len(_r) else _r
        _t = _wl
        if len(_r):
            _pf = _r['p_fdr'].iloc[0]
            _st = '***' if _pf < 0.001 else '**' if _pf < 0.01 else '*' if _pf < 0.05 else 'ns'
            _t = f"{_wl}  HR {_r['hr_per_SD_lower'].iloc[0]:.2f}, pFDR={_pf:.3f} {_st}"
        _ax.set_title(_t, fontsize=11)
        _ax.set_xlim(0, 15); _ax.set_ylim(0, 1.02)
        _ax.set_xlabel('Years from diagnosis')
        if _j == 0: _ax.set_ylabel(f'Free of {_full}')
        else: _ax.set_ylabel('')
        _ax.legend(fontsize=8, loc='lower left', frameon=False)
_plt9.tight_layout()
_pd9 = None
for _nm in ['PLOT_DIR', 'PUB_DIR', 'FIG_DIR', 'PUBLICATION_DIR']:
    if _nm in globals(): _pd9 = globals()[_nm]; break
if _pd9 is None:
    _pd9 = os.path.join(OUT_DIR, 'publication_plots')
    if not os.path.isdir(_pd9): _pd9 = OUT_DIR
_fig9.savefig(os.path.join(_pd9, 'fig_survival_early_VH_CI.png'), dpi=300, bbox_inches='tight')
_plt9.close(_fig9)
print("  Saved:", os.path.join(_pd9, 'fig_survival_early_VH_CI.png'))
