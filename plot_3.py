
import numpy as np
import pandas as pd
from scipy.stats import kurtosis

# matplotlib
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams.update({'font.family': 'serif', 'font.size': 18,
    'axes.labelsize': 20,'axes.titlesize': 24, 'figure.titlesize' : 28})
matplotlib.rcParams['text.usetex'] = True

# %%
colours = [r'#fc8d62',r'#66c2a5', r'#8da0cb', r'#e78ac3', r'#ffd92f',r'#a6d854']
labels = ['CTH','LTH','KTH','LTU','TTY','AU']

# %% Load the data
prefix = r''
df = pd.read_pickle(prefix + r'Nordic_data.pkl', compression='zip')

df = df.fillna(method='backfill')
lim = 0.025

# def clean_excursions(d,lim):
#     temp = np.array(d[:])[:]
#     temp[np.where((np.diff(temp) > lim*np.std(np.diff(temp))) | (np.diff(temp) < -lim*np.std(np.diff(temp))))[0]] = np.nan
#     return temp
#
# df['KTH'] = clean_excursions(df['KTH'].values, lim)
# df['LTH'] = clean_excursions(df['LTH'].values, lim)
# df['CTH'] = clean_excursions(df['CTH'].values, lim)
# df['Tampere'] = clean_excursions(df['Tampere'].values, lim)
# df['LTU'] = clean_excursions(df['LTU'].values, lim)
# df['Aalto'] = clean_excursions(df['Aalto'].values, lim)

df = df.fillna(method='backfill')

# Round the values to ensure similar precision on all datasets
KTH = np.round(df['KTH'].values, 4)
LTH = np.round(df['LTH'].values, 4)
CTH = np.round(df['CTH'].values, 4)
Tampere = np.round(df['Tampere'].values, 4)
LTU = np.round(df['LTU'].values, 4)
Aalto = np.round(df['Aalto'].values, 4)

# %% Variance and Kurtosis -- warning: this is a very slow process
ts = np.zeros([LTU.size, 6])
ts[:,0] = CTH
ts[:,1] = LTH
ts[:,2] = KTH
ts[:,3] = LTU
ts[:,4] = Tampere
ts[:,5] = Aalto

tau = 250
var = np.zeros((tau-1,7))
Kurt = np.zeros((tau-1,7))
for i in range(1,tau):
    var[i-1,:6]  = np.nanvar(ts[i:]-ts[:-i], axis=0)
    Kurt[i-1,:6]  = kurtosis(ts[i:]-ts[:-i], nan_policy='omit') + 3


# %% plot
fig, ax = plt.subplots(2,1, figsize=(7,6));

x = np.linspace(0.02,var.shape[0]/50,var.shape[0])
for i in range(6):
    ax[0].plot(x[:],var[:,i],color = colours[i], label=labels[i], linewidth=3)

ax[0].plot(x[160:], 0.0000036*(x[160:]**1.9), '--', color = 'black', linewidth=2)
ax[0].plot(x[160:], 0.0000105*(x[160:]**1.4), '--', color = 'black', linewidth=2)

ax[0].plot(x[:], 0.0000035*(x[:]**1), ':', color = 'black', linewidth=2)
fig.text(0.64,0.64,r'Brownian motion',fontsize=20, rotation = 4);

ax[0].set_xticklabels([])
ax[0].set_yticks([0,0.00002,0.00004,0.00006,0.00008])
ax[0].set_yticklabels([0,2,4,6,8])
fig.text(0.01,0.92,r'$\times 10^{-5}$',fontsize=18);

ax[0].set_ylabel(r'$\sigma(\Delta f_\tau)^2$',labelpad=7)
ax[0].legend(loc=2,fontsize=18, ncol=3, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2)
fig.subplots_adjust(left=0.15, bottom=0.22, right=.99, top=0.99, hspace=0.03, wspace=0.18)
fig.text(0.75,0.88,r'$\tau^{1.4}$',fontsize=26);
fig.text(0.83,0.76,r'$\tau^{1.9}$',fontsize=26);

axi = fig.add_axes([.3,.3,.35,.23])

for i in range(6):
    ax[1].semilogy(x[1:],Kurt[1:,i],color = colours[i], label=labels[i], linewidth=3)
    axi.semilogy(x[1:25],Kurt[1:25,i],color = colours[i], linewidth=3)

ax[1].semilogy(x[1:],x[1:]*0 +3, ':', color = 'black', linewidth=2)
axi.semilogy(x[1:25],x[1:25]*0 +3, ':', color = 'black', linewidth=2)

axi.set_ylim([2.4,12])
axi.set_yticks([3,4,5,6,10])
axi.set_yticklabels([3,4,'',6,10])
axi.set_xticks([0.1,.2,.3,.4,.5])

ax[1].set_ylim([2.4,None])
ax[1].set_yticks([3,4,6,10,20,40])
ax[1].set_yticklabels([3,4,6,10,20,40])
ax[1].set_ylabel(r'$\kappa(\Delta f_\tau)$')
ax[1].set_xlabel(r'$\tau$ [s]')
ax[1].legend(loc=1,fontsize=18, ncol=2, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2)
fig.subplots_adjust(left=0.11, bottom=0.12, right=.99, top=0.99, hspace=0.05, wspace=0.18)
fig.text(0.01,0.96,r'\textbf{a}',fontsize=26);
fig.text(0.01,0.51,r'\textbf{b}',fontsize=26);
fig.savefig(prefix + 'Fig_2.pdf', dpi=600, transparent = True)
