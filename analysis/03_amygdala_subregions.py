"""Group nuclei into sCLR, CMA and BLA and define the baseline PD sample."""

SCLR_raw = ['amyg_lh_Anterior_amygdaloid_area_AAA','amyg_lh_Cortical_nucleus','amyg_lh_Corticoamygdaloid_transitio',
            'amyg_rh_Anterior_amygdaloid_area_AAA','amyg_rh_Cortical_nucleus','amyg_rh_Corticoamygdaloid_transitio']
BLA_raw  = ['amyg_lh_Lateral_nucleus','amyg_lh_Basal_nucleus','amyg_lh_Accessory_Basal_nucleus','amyg_lh_Paralaminar_nucleus',
            'amyg_rh_Lateral_nucleus','amyg_rh_Basal_nucleus','amyg_rh_Accessory_Basal_nucleus','amyg_rh_Paralaminar_nucleus']
CMA_raw  = ['amyg_lh_Central_nucleus','amyg_lh_Medial_nucleus','amyg_rh_Central_nucleus','amyg_rh_Medial_nucleus']
for raw_list, name in [(SCLR_raw,'z_sCLR'),(BLA_raw,'z_BLA'),(CMA_raw,'z_CMA')]:
    lh = [f'z_{c}' for c in raw_list if '_lh_' in c and f'z_{c}' in df.columns]
    rh = [f'z_{c}' for c in raw_list if '_rh_' in c and f'z_{c}' in df.columns]
    if lh and rh: df[name] = (df[lh].mean(axis=1)+df[rh].mean(axis=1))/2
if 'z_amyg_lh_Whole_amygdala' in df.columns:
    df['z_Whole_amygdala'] = (df['z_amyg_lh_Whole_amygdala']+df['z_amyg_rh_Whole_amygdala'])/2

SUBFIELD_NAMES = ['Lateral_nucleus','Basal_nucleus','Accessory_Basal_nucleus','Paralaminar_nucleus',
                  'Central_nucleus','Medial_nucleus','Cortical_nucleus',
                  'Anterior_amygdaloid_area_AAA','Corticoamygdaloid_transitio']
for sf in SUBFIELD_NAMES:
    lz,rz = f'z_amyg_lh_{sf}', f'z_amyg_rh_{sf}'
    if lz in df.columns and rz in df.columns: df[f'z_bi_{sf}'] = (df[lz]+df[rz])/2

df['visit_rank'] = df['Visit'].map(VISIT_ORDER)
bl_all = df.sort_values('visit_rank').drop_duplicates('PATNO', keep='first').copy()

SECTORS = [('z_sCLR','sCLR'),('z_CMA','CMA'),('z_BLA','BLA')]
WHOLE = ('z_Whole_amygdala','Whole amygdala')
_SUBFIELD_LABEL_OVERRIDES = {'Anterior_amygdaloid_area_AAA': 'AAA', 'Corticoamygdaloid_transitio': 'CAT'}
NUCLEI = [(f'z_bi_{sf}', _SUBFIELD_LABEL_OVERRIDES.get(sf, sf.replace('_',' ')))
          for sf in SUBFIELD_NAMES if f'z_bi_{sf}' in df.columns]
print(f"Sectors: {[l for _,l in SECTORS]} | Whole amygdala (standalone) | Nuclei (post-hoc only): {[l for _,l in NUCLEI]}")

bl_pd = bl_all[bl_all['Group']=='PD'].copy()
bl_pd['st_II']  = np.where(bl_pd['subtype'].isna(), np.nan, (bl_pd['subtype']=='II_Intermediate').astype(float))
bl_pd['st_III'] = np.where(bl_pd['subtype'].isna(), np.nan, (bl_pd['subtype']=='III_Diffuse').astype(float))
COV = ['st_II','st_III']
n_unclassified = bl_pd['subtype'].isna().sum()
print(f"PD subjects: {len(bl_pd)} | Unclassified (excluded from subtype-adjusted models): {n_unclassified}")
print(bl_pd['subtype'].value_counts(dropna=False))
