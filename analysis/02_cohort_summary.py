"""Participant flow and follow-up duration."""

def visit_count_summary(data, group_filter_col='Group', group_val=None, label=''):
    sub = data[data[group_filter_col]==group_val] if group_val is not None else data
    sub = sub[sub['Visit'].isin(VISITS)]
    vpp = sub.groupby('PATNO')['Visit'].nunique()
    counts = {f'n_visits_ge_{k}': int((vpp >= k).sum()) for k in [1,2,3,4,5]}
    print(f"  {label}: >=1: {counts['n_visits_ge_1']} | >=2: {counts['n_visits_ge_2']} | "
          f">=3: {counts['n_visits_ge_3']} | >=4: {counts['n_visits_ge_4']} | all 5: {counts['n_visits_ge_5']}")
    return {'label':label, **counts}

flow_rows = []
print("Before QC:")
flow_rows.append(visit_count_summary(harm, 'Group', 'Control', 'HC, before QC'))
flow_rows.append(visit_count_summary(harm, 'Group', 'PD', 'PD, before QC'))
for st in SUBTYPES:
    harm_st = harm[harm['subtype']==st] if 'subtype' in harm.columns else harm.iloc[0:0]
    flow_rows.append(visit_count_summary(harm_st, None, None, f'PD Subtype {PL[st]}, before QC'))
print("After QC:")
flow_rows.append(visit_count_summary(df, 'Group', 'Control', 'HC, after QC'))
flow_rows.append(visit_count_summary(df, 'Group', 'PD', 'PD, after QC'))
for st in SUBTYPES:
    df_st = df[df['subtype']==st] if 'subtype' in df.columns else df.iloc[0:0]
    flow_rows.append(visit_count_summary(df_st, None, None, f'PD Subtype {PL[st]}, after QC'))
pd.DataFrame(flow_rows).to_csv(f'{OUT_DIR}/participant_flow.csv', index=False)

def followup_summary(data, group_label):
    sub = data[data['Visit'].isin(VISITS)].copy()
    visits_per_subj = sub.groupby('PATNO')['Visit'].nunique()
    multi_patnos = visits_per_subj[visits_per_subj >= 2].index
    n_no_followup = (visits_per_subj == 1).sum()
    n_total = visits_per_subj.shape[0]
    last_visit_years = sub[sub['PATNO'].isin(multi_patnos)].groupby('PATNO')['years_from_first_mri'].max()
    if len(last_visit_years) == 0:
        print(f"  {group_label:<20}: n_total={n_total}, n with >=2 visits=0, no follow-up data available")
        return {'group':group_label,'n_total':n_total,'n_with_followup':0,'n_no_followup':n_no_followup,
               'mean_years':np.nan,'sd_years':np.nan,'median_years':np.nan,'min_years':np.nan,'max_years':np.nan}
    mean_y = last_visit_years.mean(); sd_y = last_visit_years.std()
    median_y = last_visit_years.median(); min_y = last_visit_years.min(); max_y = last_visit_years.max()
    print(f"  {group_label:<20}: n_total={n_total}, n with >=2 visits={len(last_visit_years)} ({n_no_followup} with only baseline), "
          f"follow-up = {mean_y:.2f}+/-{sd_y:.2f} yrs (median={median_y:.2f}, range=[{min_y:.2f}, {max_y:.2f}])")
    return {'group':group_label,'n_total':n_total,'n_with_followup':len(last_visit_years),'n_no_followup':n_no_followup,
           'mean_years':mean_y,'sd_years':sd_y,'median_years':median_y,'min_years':min_y,'max_years':max_y}

followup_rows = []
followup_rows.append(followup_summary(df[df['Group']=='Control'], 'HC'))
followup_rows.append(followup_summary(df[df['Group']=='PD'], 'All PD'))
_df_subtype_lookup = df.copy()
_df_subtype_lookup['subtype'] = _df_subtype_lookup['PATNO'].map(patno_berg)
for st in SUBTYPES:
    followup_rows.append(followup_summary(_df_subtype_lookup[(_df_subtype_lookup['Group']=='PD') & (_df_subtype_lookup['subtype']==st)],
                                          f'PD Subtype {PL[st]}'))
pd.DataFrame(followup_rows).to_csv(f'{OUT_DIR}/followup_duration_summary.csv', index=False)
print(f"Saved: {OUT_DIR}/followup_duration_summary.csv")
