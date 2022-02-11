import numpy as np
import pandas as pd
# from scipy.stats import kurtosis, linregress, pearsonr

# Maps
# import geopandas
# import geoplot

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

# %%
def histograms_increment(bins, tau, lim = (0.04,0.04), *args):

    histo = np.zeros([bins, len(args), tau])
    edges = np.zeros([bins+1, len(args), tau])
    i = 0
    for el in args:
        for j in range(1,tau):
            histo[:,i,j-1], edges[:,i,j-1]  = np.histogram(el[j:]-el[:-1*j], bins=bins, range = lim, density=True)
        i +=1

    return histo, edges

histo, edges = histograms_increment(41, 51, (-0.02,0.02), CTH, LTH, KTH, LTU, Tampere, Aalto)

# %% Maps
# world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))
# Norway = world[world['name'] == "Norway"]
# Sweden = world[world['name'] == "Sweden"]
# Finland = world[world['name'] == "Finland"]
#
# lat =    [57.689769,55.711850,59.350029,65.617792,61.494200,60.186463]
# lon =    [11.973701,13.210120,18.070009,22.135986,23.780750,24.829515]
# locations_of_measures = [[lat[i],lon[i]] for i in range(6)]
# where = [(-29,-7),(0,10),(0,-22),(0,10),(0,10),(0,-22)]
# where = [(0,12),(0,-24),(0,12),(0,12),(0,10),(0,-24)]
# labels_b = [r'\textbf{CTH}',r'\textbf{LTH}',r'\textbf{KTH}',r'\textbf{LTU}',r'\textbf{TTY}',r'\textbf{AU}']

# %%
fig, ax = plt.subplots(2,3, figsize=(15,6));

for i in range(6):
    ax[1,i//2].semilogy(edges[1:,i,0],histo[:,i,0], color = colours[i], label=labels[i], linewidth=3)
    ax[1,i//2].semilogy(edges[1:,i,49],300*histo[:,i,49], '--', color = colours[i], linewidth=3)

var = np.var(CTH[1:]-CTH[:-1])
ax[1,0].semilogy((edges[1:,0,0] + edges[:-1,0,0])/2,np.exp(-(edges[:-1,0,0]**2)/(2*var))/(np.sqrt(var*2*np.pi)), ':', color = 'black', linewidth=2)
var = np.var(KTH[1:]-KTH[:-1])
ax[1,1].semilogy((edges[1:,2,0] + edges[:-1,2,0])/2,np.exp(-(edges[:-1,2,0]**2)/(2*var))/(np.sqrt(var*2*np.pi)), ':', color = 'black', linewidth=2)
var = np.var(Tampere[1:]-Tampere[:-1])
ax[1,2].semilogy((edges[1:,4,0] + edges[:-1,4,0])/2,np.exp(-(edges[:-1,4,0]**2)/(2*var))/(np.sqrt(var*2*np.pi)), ':', color = 'black', linewidth=2)

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], linestyle='-', color='black',  lw=1.5, label=r'$\tau=0.02\!~$s'),
                Line2D([0], [0], linestyle='--', color='black',  lw=1.5, label=r'$\tau=1\!~$s')]

lab1 = [r'$\tau=0.02\!~$s', r'$\tau=1\!~$s']

leg1 = ax[1,0].legend(custom_lines, lab1, loc=4, handlelength=0.7, handletextpad=0.3, ncol = 2, columnspacing=0.5)
leg2 = ax[1,1].legend(custom_lines, lab1, loc=4, handlelength=0.7, handletextpad=0.3, ncol = 2, columnspacing=0.5)
leg3 = ax[1,2].legend(custom_lines, lab1, loc=4, handlelength=0.7, handletextpad=0.3, ncol = 2, columnspacing=0.5)

ax[1,0].add_artist(leg1)
ax[1,1].add_artist(leg2)
ax[1,2].add_artist(leg3)

locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(.2,.4,.6,.8),numticks=12)
ax[1,0].yaxis.set_minor_locator(locmin)
ax[1,0].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[1,1].yaxis.set_minor_locator(locmin)
ax[1,1].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[1,2].yaxis.set_minor_locator(locmin)
ax[1,2].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

# world.plot(ax = ax[0,0], color='#D8D8D8',  edgecolor='white');
# Norway.plot(ax = ax[0,0], color='#A0A0A0',  edgecolor='white');
# Sweden.plot(ax = ax[0,0], color='#A0A0A0',  edgecolor='white');
# Finland.plot(ax = ax[0,0], color='#A0A0A0',  edgecolor='white');

# for i in range(6):
#     ax[0,0].scatter(locations_of_measures[i][1],locations_of_measures[i][0], color = colours[i], s=150, edgecolors = 'black')
#     ax[0,0].annotate(labels_b[i], (locations_of_measures[i][1],locations_of_measures[i][0]),textcoords="offset points", xytext=where[i], ha='center',fontsize=18)

ax[0,1].plot(    CTH[:1800*50]*1000 + 50*4,color = colours[0], label=labels[0], linewidth=3)
ax[0,1].plot(    LTH[:1800*50]*1000 + 30*4,color = colours[1], label=labels[1], linewidth=3)
ax[0,1].plot(    KTH[:1800*50]*1000 + 10*4,color = colours[2], label=labels[2], linewidth=3)
ax[0,1].plot(    LTU[:1800*50]*1000 - 10*4,color = colours[3], label=labels[3], linewidth=3)
ax[0,1].plot(Tampere[:1800*50]*1000 - 30*4,color = colours[4], label=labels[4], linewidth=3)
ax[0,1].plot(  Aalto[:1800*50]*1000 - 50*4,color = colours[5], label=labels[5], linewidth=3)

ax[0,2].plot(    CTH[1350*50:1350*50 + 100]*1000 + 25*1,color = colours[0], label=labels[0], linewidth=3)
ax[0,2].plot(    LTH[1350*50:1350*50 + 100]*1000 + 15*1,color = colours[1], label=labels[1], linewidth=3)
ax[0,2].plot(    KTH[1350*50:1350*50 + 100]*1000 + 5*1,color = colours[2], label=labels[2], linewidth=3)
ax[0,2].plot(    LTU[1350*50:1350*50 + 100]*1000 - 5*1,color = colours[3], label=labels[3], linewidth=3)
ax[0,2].plot(Tampere[1350*50:1350*50 + 100]*1000 - 15*1,color = colours[4], label=labels[4], linewidth=3)
ax[0,2].plot(  Aalto[1350*50:1350*50 + 100]*1000 - 25*1,color = colours[5], label=labels[5], linewidth=3)

ax[0,1].set_ylim([-380,320])
ax[0,1].set_xticks([0,450*50,900*50,1350*50,1800*50])
ax[0,1].set_xticklabels([0,450,900,1350,1800])
ax[0,2].set_ylim([-24,44])
ax[0,2].set_yticks([-15,0,15,30])
ax[0,2].set_xticks([0,20,40,60,80,100])
ax[0,2].set_xticklabels([0,20*20,40*20,60*20,80*20,100*20])

ax[0,1].set_ylabel(r'$f$ [mHz]',labelpad=-14)
ax[0,1].set_xlabel(r'$t$ [s]')
ax[0,2].set_xlabel(r'$t$ [ms]')
ax[1,0].set_ylabel(r'PDF',labelpad=-8)

ax[0,1].annotate('', xy=(1350*50, 220),xytext=(1350*50, 310),
    horizontalalignment="center", arrowprops=dict(arrowstyle='->',lw=1),fontsize=26)
ax[0,1].annotate('', xy=(1350*50, -210),xytext=(1350*50, -300),
    horizontalalignment="center", arrowprops=dict(arrowstyle='->',lw=1),fontsize=26)

# ax[0,0].set_ylabel(r'Latitude')
# ax[0,0].set_xlabel(r'Longitude')
# ax[0,0].set_xlim([-1,37]); ax[0,0].set_ylim([52,71.5]);

[ax[1,i].set_xlabel(r'$\Delta f_\tau$') for i in range(3)]

fig.text(0.002,0.95,r'\textbf{a}',fontsize=28);
fig.text(0.335,0.95,r'\textbf{b}',fontsize=28);
fig.text(0.665,0.95,r'\textbf{c}',fontsize=28);
fig.text(0.002,0.45,r'\textbf{d}',fontsize=28);
fig.text(0.335,0.45,r'\textbf{e}',fontsize=28);
fig.text(0.665,0.45,r'\textbf{f}',fontsize=28);

ax[1,0].set_ylim([2e-4,2e5])
ax[1,1].set_ylim([6e-5,2e5])
ax[1,2].set_ylim([2e-4,2e5])

[ax[1,i].legend(loc=2,fontsize=17, ncol=1, handlelength=0.5, columnspacing=0.5, handletextpad = 0.2) for i in range(3)]
fig.subplots_adjust(left=0.05, bottom=0.11, right=.99, top=0.99, hspace=0.32, wspace=0.19)

# fig.savefig(prefix + 'Fig_1_no_map.pdf', dpi=600, transparent = True)
