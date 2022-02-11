import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d

# matplotlib
import matplotlib
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams.update({'font.family': 'serif', 'font.size': 18,
    'axes.labelsize': 20,'axes.titlesize': 24, 'figure.titlesize' : 28})
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

colours = [r'#fc8d62',r'#66c2a5', r'#8da0cb', r'#e78ac3', r'#ffd92f',r'#a6d854']
labels = ['CTH','LTH','KTH','LTU','TTY','AU']

prefix = r''
df = pd.read_pickle(prefix + r'Nordic_data.pkl', compression='zip')

# Remove the few NaNs left
df = df.fillna(method='backfill')
lim = 0.025

def clean_excursions(d,lim):
    d[:-1][np.diff(d)>lim] = np.nan
    d[:-1][np.diff(d)<-1*lim] = np.nan
    return d

df['KTH'] = clean_excursions(df['KTH'].values, lim)
df['LTH'] = clean_excursions(df['LTH'].values, lim)
df['CTH'] = clean_excursions(df['CTH'].values, lim)
df['Tampere'] = clean_excursions(df['Tampere'].values, lim)
df['LTU'] = clean_excursions(df['LTU'].values, lim)
df['Aalto'] = clean_excursions(df['Aalto'].values, lim)

df = df.fillna(method='backfill')

# Round the values to ensure similar precision on all datasets
KTH = np.round(df['KTH'].values, 3)
LTH = np.round(df['LTH'].values, 3)
CTH = np.round(df['CTH'].values, 3)
Tampere = np.round(df['Tampere'].values, 3)
LTU = np.round(df['LTU'].values, 3)
Aalto = np.round(df['Aalto'].values, 3)

# %% Pearson - Warning: very slow process, data is included for expidiency
# ts = np.zeros([LTU.size, 6])
# ts[:,0] = CTH; ts[:,1] = LTH; ts[:,2] = KTH; ts[:,3] = LTU; ts[:,4] = Tampere; ts[:,5] = Aalto;
#
# pear = np.zeros((6,6,1000))
#
# # %%
# for k in tqdm(range(1,1001)):
#     increments = (np.asarray(ts[k:]) - np.asarray(ts[:-k]))[:]*1000
#     for i in range(6):
#         for j in range(6):
#             if i!=j:
#                 pear[i,j,k-1] = pearsonr(increments[:,i],increments[:,j])[0]

## this is the data for the commented script above
pear = np.load(prefix + 'pearson.npz')['pear']

# %%
x0, k, b = np.zeros((6,6)), np.zeros((6,6)), np.zeros((6,6))
x_ = np.linspace(0.02,(pear.shape[2]-2)/50,pear.shape[2]-2)

pear_clean = np.zeros_like(pear)

for i in range(6):
    for j in range(6):
        if i!=j:
            pear_clean[i,j,:] = gaussian_filter1d(pear[i,j,:], sigma=3)

# %%
distances = np.zeros((6,6))

### Using openstreetmaps and OSRM: http://map.project-osrm.org/
# CTH to
distances[2,5] = 2141.2; distances[5,2] = 2141.2;
distances[2,3] = 2018.6; distances[3,2] = 2018.6;
distances[2,4] = 1255.2; distances[4,2] = 1255.2;
distances[2,1] = 263.3;  distances[1,2] = 263.3;
distances[2,0] = 436.4;  distances[0,2] = 436.4;

# KTH to
distances[0,1] = 569.1;  distances[1,0] = 569.1;
distances[0,3] = 1679.6; distances[3,0] = 1679.6;
distances[0,4] = 932.0;  distances[4,0] = 932.0;
distances[0,5] = 1802.1; distances[5,0] = 1802.1;

# LTU to
distances[4,1] = 1461.0; distances[1,4] = 1461.0;
distances[4,3] = 754.1;  distances[3,4] = 754.1;
distances[4,5] = 876.5;  distances[5,4] = 876.5;

# LTH to
distances[1,3] = 2248.3; distances[3,1] = 2248.3;
distances[1,5] = 2370.8;  distances[5,1] = 2370.8;

# Aalto to
distances[5,3] = 180.5; distances[3,5] = 180.5;


# %%
pear.shape
lin_dist = np.triu(distances).flatten()
lin_dist = lin_dist[np.abs(lin_dist)!=0.]

lin_pear = np.zeros((15,1000))
for i in range(1,1000):
    _ = np.triu(pear_clean[...,i-1]).flatten()
    lin_pear[:,i-1] = _[np.abs(_)!=0.]

# %% fits
def fun3(x, a,c):
    return a*x**c

sorting = np.argsort(lin_dist[:-3])
local_lin_dist = lin_dist[sorting]
local_lin_pear = lin_pear[sorting,5]

positives = local_lin_pear>0
fit_short3, error_short3 = curve_fit(fun3, local_lin_dist[positives], np.log(local_lin_pear[positives]))

# %%
arr = np.linspace(0.02,20,1000)
dis = np.zeros(15)
for i in range(15):
    dis[i] = arr[:-1][np.diff(lin_pear[i,:]>0.5)][0]

def fun2(x, a, b):
    return a*x + b

sorting = np.argsort(lin_dist[:-3])

fit_medium, error_medium = curve_fit(fun2, lin_dist[sorting], dis[sorting])
fit_medium_std = np.std(dis[sorting] - fun2(lin_dist[sorting], *fit_medium))

# %% correlations in space plot
fig, ax = plt.subplots(2,1, figsize=(7,6));

gradient = matplotlib.cm.get_cmap('cividis', 256)

x = np.linspace(200,2400,220)
ax[0].plot(x, np.exp(fun3(x, *fit_short3)), '-.', color='black', lw=2,label=r'Fit $\tau=0.1~\!$s: exp($a x^k$)')
x_ = np.linspace(100,2500,220)
ax[0].plot(x_, x_*0, ':', color='black', lw=1)

ax[0].plot(lin_dist[[0,1,2,3,4,5,6,7,8,9,10,11]], lin_pear[[0,1,2,3,4,5,6,7,8,9,10,11],5], 'o', ms=10, markeredgecolor = 'black', color=gradient(((i+16)/16)))
ax[0].plot(lin_dist[-3:], lin_pear[-3:,5], 'o', ms=10, markeredgecolor = gradient(((1)/16)), markerfacecolor='none', label=r'Excluded from fit')

x = np.linspace(200,2400,220)
ax[1].plot(x, fun2(x, *fit_medium), '-.', color='black', lw=2,label=r'Fit $a x + b$')
ax[1].plot(lin_dist[:-3], dis[:-3], 'D', ms=8, markeredgecolor = 'black', color=gradient(((3)/16)))
ax[1].plot(lin_dist[-3:], dis[-3:], 'D', ms=8, markeredgecolor = 'black', color=gradient(((3)/16)), markerfacecolor='none', label=r'Excluded from fit')

ax[1].fill_between(x, fun2(x, *fit_medium), fun2(x, *fit_medium)+fit_medium_std, color='black', alpha=0.1)
ax[1].fill_between(x, fun2(x, *fit_medium), fun2(x, *fit_medium)-fit_medium_std, color='black', alpha=0.1)
x_ = np.linspace(100,2500,220)

ax[0].set_xticks([0,400,800,1200,1600,2000,2400])
ax[0].set_xticklabels([])
ax[0].set_xlim([-30,2550])
ax[0].set_ylim([-0.08,.45])
ax[0].set_ylabel(r'$c(\Delta f_\tau^X,\Delta f_\tau^{Y})$')

ax[1].set_xlim([-30,2550])
ax[1].set_xticks([0,400,800,1200,1600,2000,2400])
ax[1].set_ylim([-0.08,2.1])
ax[1].set_yticks([0,0.5,1,1.5,2])
ax[1].set_ylabel(r'$\tau_{c>0.5}$')
ax[1].set_xlabel(r'Distance [km]')

fig.text(0.02,0.95,r'\textbf{a}',fontsize=26);
fig.text(0.02,0.51,r'\textbf{b}',fontsize=26);
ax[0].legend(loc=1,fontsize=18, ncol=1, handlelength=1.15, columnspacing=0.5, handletextpad = 0.2)
ax[1].legend(loc=2,fontsize=18, ncol=1, handlelength=1.15, columnspacing=0.5, handletextpad = 0.2)
fig.subplots_adjust(left=0.15, bottom=0.12, right=.99, top=0.99, hspace=0.05, wspace=0.18)
fig.savefig(prefix + 'Fig_5.pdf', dpi=600, transparent = True)
