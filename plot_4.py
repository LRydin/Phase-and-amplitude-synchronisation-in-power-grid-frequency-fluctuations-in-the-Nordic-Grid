import numpy as np
import pandas as pd
from scipy.stats import kurtosis, pearsonr
from scipy.optimize import curve_fit

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

ts = np.zeros([LTU.size, 6])
ts[:,0] = CTH; ts[:,1] = LTH; ts[:,2] = KTH; ts[:,3] = LTU; ts[:,4] = Tampere; ts[:,5] = Aalto;

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

def sigmoid(x, x0, k, b):
    return 1 / (1 + np.exp(-k*(x-x0)))+b

for i in range(6):
    for j in range(6):
        if i!=j:
            x0[i,j], k[i,j], b[i,j] = curve_fit(sigmoid, x_, pear[i,j,:-2], method='dogbox')[0]

# %%
fig = plt.figure(figsize=(7,5))
gs = fig.add_gridspec(6, 3)

ax1 = fig.add_subplot(gs[:2, 0])
ax2 = fig.add_subplot(gs[:2, 1])
ax3 = fig.add_subplot(gs[:2, 2])
ax4 = fig.add_subplot(gs[3:, :])

increments1 = (np.asarray(ts[1:]) - np.asarray(ts[:-1]))[:]*1000
increments2 = (np.asarray(ts[20:]) - np.asarray(ts[:-20]))[:]*1000
increments3 = (np.asarray(ts[100:]) - np.asarray(ts[:-100]))[:]*1000

ax1.set_ylabel(r'$\Delta f_\tau^{\mathrm{CTH}}$',labelpad=-5)
p1 = ax1.scatter(increments1[:,0],increments1[:,5], s = 2 ,color = colours[5], alpha=0.01)
p2 = ax2.scatter(increments2[:,0],increments2[:,5], s = 2 ,color = colours[5], alpha=0.01)
p3 = ax3.scatter(increments3[:,0],increments3[:,5], s = 2 ,color = colours[5], alpha=0.01)
p1.set_rasterized(True)
p2.set_rasterized(True)
p3.set_rasterized(True)

fig.text(0.157+0.04,0.945,r'$\tau\!\!=\!0.02\!~$s',fontsize=16,ha='left')
fig.text(0.463+0.04,0.945,r'$\tau\!\!=\!0.40\!~$s',fontsize=16,ha='left')
fig.text(0.765+0.04,0.945,r'$\tau\!\!=\!2.00\!~$s',fontsize=16,ha='left')

ax1.set_xlabel(r'$\Delta f_\tau^{\mathrm{AU}}$',labelpad=4)
ax2.set_xlabel(r'$\Delta f_\tau^{\mathrm{AU}}$',labelpad=4)
ax3.set_xlabel(r'$\Delta f_\tau^{\mathrm{AU}}$',labelpad=4)
ax1.set_xlim([-22,22]); ax1.set_ylim([-22,22]);
ax2.set_xlim([-22,22]); ax2.set_ylim([-22,22]);
ax3.set_xlim([-22,22]); ax3.set_ylim([-22,22]);

ax2.yaxis.set_visible(False);
ax3.yaxis.set_visible(False);

for i in [1,2,3,4,5]:
    ax4.semilogx(x_,pear[0,i,:-2],color = colours[i], label=labels[i], linewidth=3)

ax4.semilogx([0.02],pear[0,5,0], 'o', markeredgecolor='black',markerfacecolor=colours[5],ms=10)
ax4.semilogx([0.4],pear[0,5,19], 's', markeredgecolor='black',markerfacecolor=colours[5],ms=10)
ax4.semilogx([2.00],pear[0,5,99], 'D', markeredgecolor='black',markerfacecolor=colours[5],ms=10)

ax1.plot([-18],[17], 'o', markeredgecolor='black',markerfacecolor=colours[5],ms=10)
ax2.plot([-18],[17], 's', markeredgecolor='black',markerfacecolor=colours[5],ms=10)
ax3.plot([-18],[17], 'D', markeredgecolor='black',markerfacecolor=colours[5],ms=10)

fig.text(0.3,0.48,r'\textbf{CTH \textit{vs}}',fontsize=28,ha='center',color=colours[0])

ax4.set_ylabel(r'$c(\Delta f_\tau^X,\Delta f_\tau^{\mathrm{CTH}})$')
ax4.set_xlabel(r'$\tau$ [s]')
ax4.legend(loc=4,fontsize=17, ncol=1, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2)
fig.subplots_adjust(left=0.15, bottom=0.14, right=.99, top=0.99, hspace=0, wspace=0.3)
fig.text(0.05,0.95,r'\textbf{a}',fontsize=26);
fig.text(0.39,0.95,r'\textbf{b}',fontsize=26);
fig.text(0.695,0.95,r'\textbf{c}',fontsize=26);
fig.text(0.05,0.59,r'\textbf{d}',fontsize=26);
# fig.savefig(prefix + 'Fig_4.pdf', dpi=600, transparent = True)


# %% ############################### extra plots  ##############################

# %% CTH vs
fig, ax = plt.subplots(1,1, figsize=(7,3));

for i in [1,2,3,4,5]:
    ax.semilogx(x_,pear[0,i,:-2],color = colours[i], label=labels[i], linewidth=3)

fig.text(0.3,0.83,r'\textbf{CTH \textit{vs}}',fontsize=28,ha='center',color=colours[0])

ax.set_ylabel(r'$c(\Delta f_\tau,\Delta f_\tau^{\mathrm{CTH}})$')
ax.set_xlabel(r'$\tau$ [s]')
ax.legend(loc=4,fontsize=17, ncol=1, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2)
fig.subplots_adjust(left=0.15, bottom=0.22, right=.99, top=0.99, hspace=0.03, wspace=0.18)
# fig.savefig(prefix + 'Fig_CTH.pdf', dpi=600, transparent = True)

# %% LTH vs
fig, ax = plt.subplots(1,1, figsize=(7,3));

for i in [0,2,3,4,5]:
    ax.semilogx(x_,pear[1,i,:-2],color = colours[i], label=labels[i], linewidth=3)

fig.text(0.3,0.83,r'\textbf{LTH \textit{vs}}',fontsize=28,ha='center',color=colours[1])

ax.set_ylabel(r'$c(\Delta f_\tau,\Delta f_\tau^{\mathrm{LTH}})$')
ax.set_xlabel(r'$\tau$ [s]')
ax.legend(loc=4,fontsize=17, ncol=1, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2)
fig.subplots_adjust(left=0.15, bottom=0.22, right=.99, top=0.99, hspace=0.03, wspace=0.18)
# fig.savefig(prefix + 'Fig_LTH.pdf', dpi=600, transparent = True)

# %% KTH vs
fig, ax = plt.subplots(1,1, figsize=(7,3));

for i in [0,1,3,4,5]:
    ax.semilogx(x_,pear[2,i,:-2],color = colours[i], label=labels[i], linewidth=3)

fig.text(0.3,0.83,r'\textbf{KTH \textit{vs}}',fontsize=28,ha='center',color=colours[2])

ax.set_ylabel(r'$c(\Delta f_\tau,\Delta f_\tau^{\mathrm{KTH}})$')
ax.set_xlabel(r'$\tau$ [s]')
ax.legend(loc=4,fontsize=17, ncol=1, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2)
fig.subplots_adjust(left=0.15, bottom=0.22, right=.99, top=0.99, hspace=0.03, wspace=0.18)
# fig.savefig(prefix + 'Fig_KTH.pdf', dpi=600, transparent = True)

# %% LTU vs
fig, ax = plt.subplots(1,1, figsize=(7,3));

for i in [0,1,2,4,5]:
    ax.semilogx(x_,pear[3,i,:-2],color = colours[i], label=labels[i], linewidth=3)

fig.text(0.3,0.83,r'\textbf{LTU \textit{vs}}',fontsize=28,ha='center',color=colours[3])

ax.set_ylabel(r'$c(\Delta f_\tau,\Delta f_\tau^{\mathrm{LTU}})$')
ax.set_xlabel(r'$\tau$ [s]')
ax.legend(loc=4,fontsize=17, ncol=1, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2)
fig.subplots_adjust(left=0.15, bottom=0.22, right=.99, top=0.99, hspace=0.03, wspace=0.18)
# fig.savefig(prefix + 'Fig_LTU.pdf', dpi=600, transparent = True)

# %% Tampere vs
fig, ax = plt.subplots(1,1, figsize=(7,3));

for i in [0,1,2,3,5]:
    ax.semilogx(x_,pear[4,i,:-2],color = colours[i], label=labels[i], linewidth=3)

fig.text(0.35,0.83,r'\textbf{Tampere \textit{vs}}',fontsize=28,ha='center',color=colours[4])

ax.set_ylabel(r'$c(\Delta f_\tau,\Delta f_\tau^{\mathrm{Tampere}})$')
ax.set_xlabel(r'$\tau$ [s]')
ax.legend(loc=4,fontsize=17, ncol=1, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2)
fig.subplots_adjust(left=0.15, bottom=0.22, right=.99, top=0.99, hspace=0.03, wspace=0.18)
# fig.savefig(prefix + 'Fig_TTY.pdf', dpi=600, transparent = True)

# %% Aalto vs
fig, ax = plt.subplots(1,1, figsize=(7,3));

for i in [0,1,2,3,4]:
    ax.semilogx(x_,pear[5,i,:-2],color = colours[i], label=labels[i], linewidth=3)

fig.text(0.3,0.83,r'\textbf{Aalto \textit{vs}}',fontsize=28,ha='center',color=colours[5])

ax.set_ylabel(r'$c(\Delta f_\tau,\Delta f_\tau^{\mathrm{Aalto}})$')
ax.set_xlabel(r'$\tau$ [s]')
ax.legend(loc=4,fontsize=17, ncol=1, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2)
fig.subplots_adjust(left=0.15, bottom=0.22, right=.99, top=0.99, hspace=0.03, wspace=0.18)
# fig.savefig(prefix + 'Figs/pearson/Fig_AU.pdf', dpi=600, transparent = True)

################################# full figure ##################################

# %% CTH vs
fig, ax = plt.subplots(2,3, figsize=(15,2.8125*2));

for i in [1,2,3,4,5]:
    ax[0,0].semilogx(x_,pear[0,i,:-2],color = colours[i], label=labels[i], linewidth=3)

ax[0,0].text(0.02,0.8,r'\textbf{CTH \textit{vs}}',fontsize=28,ha='left',color=colours[0])

for i in [0,2,3,4,5]:
    ax[0,1].semilogx(x_,pear[1,i,:-2],color = colours[i], label=labels[i], linewidth=3)

ax[0,1].text(0.02,0.8,r'\textbf{LTH \textit{vs}}',fontsize=28,ha='left',color=colours[1])

for i in [0,1,3,4,5]:
    ax[0,2].semilogx(x_,pear[2,i,:-2],color = colours[i], label=labels[i], linewidth=3)

ax[0,2].text(0.02,0.8,r'\textbf{KTH \textit{vs}}',fontsize=28,ha='left',color=colours[2])

for i in [0,1,2,4,5]:
    ax[1,0].semilogx(x_,pear[3,i,:-2],color = colours[i], label=labels[i], linewidth=3)

ax[1,0].text(0.02,0.8,r'\textbf{LTU \textit{vs}}',fontsize=28,ha='left',color=colours[3])

for i in [0,1,2,3,5]:
    ax[1,1].semilogx(x_,pear[4,i,:-2],color = colours[i], label=labels[i], linewidth=3)

ax[1,1].text(0.02,0.8,r'\textbf{TTY$\!~$\textit{vs}}',fontsize=28,ha='left',color=colours[4])

for i in [0,1,2,3,4]:
    ax[1,2].semilogx(x_,pear[5,i,:-2],color = colours[i], label=labels[i], linewidth=3)

ax[1,2].text(0.02,0.8,r'\textbf{AU \textit{vs}}',fontsize=28,ha='left',color=colours[5])

[ax[0,j].set_xticks([]) for j in range(3)]

[ax[i,0].set_ylabel(r'$c(\Delta f_\tau^{X},\Delta f_\tau^{Y})$') for i in range(2)]
[ax[1,j].set_xlabel(r'$\tau$ [s]') for j in range(3)]
[[ax[i,j].legend(loc=4,fontsize=17, ncol=1, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2) for i in range(2)] for j in range(3)]
fig.subplots_adjust(left=0.07, bottom=0.12, right=.99, top=0.99, hspace=0.05, wspace=0.18)
fig.text(0.01,0.95,r'\textbf{a}',fontsize=26); fig.text(0.35,0.95,r'\textbf{b}',fontsize=26); fig.text(0.675,0.95,r'\textbf{c}',fontsize=26);
fig.text(0.01,0.50,r'\textbf{d}',fontsize=26); fig.text(0.35,0.50,r'\textbf{e}',fontsize=26); fig.text(0.675,0.50,r'\textbf{f}',fontsize=26);
# fig.savefig(prefix + 'Fig_all.pdf', dpi=600, transparent = True)
