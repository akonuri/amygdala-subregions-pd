import pandas as pd
import numpy as np
import os
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from scipy.spatial.distance import mahalanobis
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from scikit_posthocs import posthoc_conover
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

def safe_read_csv(path, retries=3, **kwargs):
    """Read a CSV via an in-memory buffer (avoids file-handle errors on network drives)."""
    import io
    last_err = None
    for attempt in range(retries):
        try:
            with open(path, 'rb') as f:
                raw = f.read()
            return pd.read_csv(io.BytesIO(raw), **kwargs)
        except OSError as e:
            last_err = e
            if attempt < retries - 1:
                print(f"  [safe_read_csv] Retry {attempt+1}/{retries} after OSError on {path}: {e}")
                continue
            raise last_err

BASE_DIR = os.environ.get('AMYG_PD_DATA_DIR')
if not BASE_DIR:
    raise SystemExit('Set AMYG_PD_DATA_DIR to the folder holding the input files (see README).')
OUT_DIR = os.environ.get('AMYG_PD_OUTPUT_DIR', os.path.join(BASE_DIR, 'amygdala_subtypes_manuscript'))
os.makedirs(OUT_DIR, exist_ok=True)
PLOT_DIR = f'{OUT_DIR}/publication_plots'
os.makedirs(PLOT_DIR, exist_ok=True)
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 11, 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.edgecolor': '#444444', 'axes.linewidth': 0.9,
})

VISITS = ['BL','V04','V06','V08','V10']
VISIT_ORDER = {'BL':0,'V04':1,'V06':2,'V08':3,'V10':4}
SUBTYPES = ['I_Mild','II_Intermediate','III_Diffuse']
GROUPS = ['HC','I_Mild','II_Intermediate','III_Diffuse']
PL = {'I_Mild':'I','II_Intermediate':'II','III_Diffuse':'III'}

EVENT_MAP = {'SC':'BL','BL':'BL','ST':'BL',
             'V02':'V04','V04':'V04',
             'V05':'V06','V06':'V06',
             'V07':'V08','V08':'V08',
             'V09':'V10','V10':'V10','V11':'V10'}

def sig(p):
    return '***' if p<0.001 else('**' if p<0.01 else('*' if p<0.05 else('+' if p<0.10 else 'ns')))
def get_grp(d, grp):
    if grp=='HC': return d[d['Group']=='Control']
    elif grp=='All_PD': return d[d['Group']=='PD']
    else: return d[d['subtype']==grp]
