import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from MFDFA import MFDFA

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

# %%
colours = [r'#fc8d62',r'#66c2a5', r'#8da0cb', r'#e78ac3', r'#ffd92f',r'#a6d854']
labels = ['CTH','LTH','KTH','LTU','TTY','AU']

# Load the data
prefix = r''
df = pd.read_pickle(prefix + r'Nordic_data.pkl', compression='zip')

# Remove the few NaNs left
df = df.fillna(method='backfill')

def clean_excursions(d,lim):
    temp = np.array(d[:])[:]
    temp[np.where((np.diff(temp) > lim*np.std(np.diff(temp))) | (np.diff(temp) < -lim*np.std(np.diff(temp))))[0]] = np.nan
    return temp

df['KTH'] = clean_excursions(df['KTH'].values, 10)
df['LTH'] = clean_excursions(df['LTH'].values, 10)
df['CTH'] = clean_excursions(df['CTH'].values, 10)
df['Tampere'] = clean_excursions(df['Tampere'].values, 10)
df['LTU'] = clean_excursions(df['LTU'].values, 10)
df['Aalto'] = clean_excursions(df['Aalto'].values, 1)

# Round the values to ensure similar precision on all datasets
KTH = np.round(df['KTH'].values, 4)
LTH = np.round(df['LTH'].values, 4)
CTH = np.round(df['CTH'].values, 4)
Tampere = np.round(df['Tampere'].values, 4)
LTU = np.round(df['LTU'].values, 4)
Aalto = np.round(df['Aalto'].values, 4)

# %% MFDFA - Warning: very slow process, data is included for expidiency
# def dfa(lag, q, *args):
#     mdfa = np.zeros([lag.size, len(q), len(args)])
#     i = 0
#     for el in args:
#         mdfa[:,:,i] = MFDFA(el, lag, q=q)[1]
#         i +=1
#
#     return mdfa

# lag =  np.linspace(50,1000,951).astype(int)
# q = [2]
# mdfa = dfa(lag, q, CTH[~np.isnan(CTH)], LTH[~np.isnan(LTH)], KTH[~np.isnan(KTH)], LTU[~np.isnan(LTU)], Tampere[~np.isnan(Tampere)], Aalto[~np.isnan(Aalto)])

## this is the data for the commented script above
with np.load('MFDFA_data.npz') as data:
     mdfa = data['mdfa']
     lag = data['lag']


# %%
# Distances between the locations
# CTH, LTH, KTH, LTU, Tampere, Aalto
distances = np.zeros(5)
distances[0] = 263.3 # LTH
distances[1] = 436.4 # KTH
distances[2] = 1255.2 # LTU
distances[3] = 2018.6 # Tampere
distances[4] = 2141.2 # Aalto

eta = np.zeros((lag[:].size,6))
for i in range(6):
    eta[:,i] = (mdfa[:,0,i] - mdfa[:,0,0])/mdfa[:,0,0]

alpha = np.zeros(6)
for i in range(6):
    alpha[i] = np.polyfit(np.log(lag[200:450]), np.log(mdfa[200:450,0,i]),1)[0]

dist = np.zeros(5)
for i in range(1,6):
    dist[i-1] = lag[:][eta[:,i] < eta.max()/10][0]/50 - 1

def fun(x,a,b):
    return a*x**2 + b

fit, error = curve_fit(fun, distances, dist)

# %% Comparative MFDFA plots and distances
fig, ax = plt.subplots(3,1, figsize=(7,10));

for i in range(6):
    ax[0].loglog(lag[:]/50,1e5*mdfa[:,0,i]/(lag[:]**(1.5)), '-', color = colours[i], label=labels[i], linewidth=3, alpha=1)

ax[0].loglog(lag[:]/50,lag[:]*0 + 1e5*(mdfa[:,0,0]/(lag[:]**(1.5))).min(), ':', color = 'black', linewidth=2, alpha=1)
ax[0].loglog(lag[200:450]/50,3.6e-1*lag[200:450]**(.38), '--', color = 'black', linewidth=2, alpha=1)


ax[0].set_yticks([2,3,4,5])
ax[0].set_yticklabels([2,3,4,5])
ax[0].set_xticks([1,2,3,5,10,15])
ax[0].set_xticklabels([1,2,3,5,10,15])
ax[0].set_ylabel(r'$F(r)/F(r)^{\mathrm{Bm}}$')
ax[0].set_xlabel(r'$r$ [s]')

fig.text(0.64,0.76,r'Brownian motion',fontsize=20);
fig.text(0.64,0.84,r'$r^{\alpha^\prime}$',fontsize=28);
fig.text(0.015,0.725,r'$\times 10^{-5}$',fontsize=18);

for i in range(6):
    ax[1].semilogx(lag[:]/50, eta[:,i], '-', color = colours[i], label=labels[i], linewidth=3, alpha=1)

ax[1].semilogx(lag[:]/50,lag[:]*0+eta.max()/10, '-', color = 'gray', linewidth=2, alpha=1)
ax[1].set_ylim([None,0.27])

for i in range(5):
    ax[1].semilogx(dist[i]+1, eta.max()/10, marker='o', markeredgecolor='black', color = colours[i+1],markersize=12)

x = np.linspace(0, distances.max(),2500)

def fun_star(x,a):
    return a*x**2

ax[2].plot(x,fun(x,*fit), '--', color='black')

for i in range(5):
    ax[2].plot(distances[i],dist[i], marker='o', markeredgecolor='black', color = colours[i+1],  label=labels[i+1],markersize=12)

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], linestyle='--', color='black',  lw=1.5),
                Line2D([0], [0], linestyle='-.', color='black',  lw=1.5)]

lab = [r'fit: $a x^{2} + b$']
leg1 = ax[2].legend(custom_lines, lab, loc=2, handlelength=0.85,
                    handletextpad=0.5, ncol = 1, columnspacing=0.5,
                    bbox_to_anchor=(0,0.75,0,0))

ax[2].add_artist(leg1)

ax[1].set_xticks([1,2,3,5,10,15])
ax[1].set_xticklabels([1,2,3,5,10,15])
ax[1].set_ylabel(r'$\eta(r)^{\mathrm{CTH}}$')
ax[1].set_xlabel(r'$r$ [s]')
ax[2].set_ylabel(r'$\chi^{\mathrm{CTH}}$')
ax[2].set_xlabel(r'Distance [km]')

ax[0].legend(loc=2,fontsize=18, ncol=3, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2)
ax[1].legend(loc=1,fontsize=18, ncol=3, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2)
ax[2].legend(loc=2,fontsize=18, ncol=5, handlelength=0.2, columnspacing=0.5, handletextpad = 0.6, markerscale=0.7)
fig.subplots_adjust(left=0.12, bottom=0.07, right=.99, top=0.99, hspace=0.30, wspace=0.18)
fig.text(0.01,0.98,r'\textbf{a}',fontsize=26);
fig.text(0.01,0.63,r'\textbf{b}',fontsize=26);
fig.text(0.01,0.31,r'\textbf{c}',fontsize=26);
fig.savefig(prefix + 'Fig_6.pdf', dpi=600, transparent = True)
