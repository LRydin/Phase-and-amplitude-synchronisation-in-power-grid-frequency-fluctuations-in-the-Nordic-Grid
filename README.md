# Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid
Related python code to the scientific publication "Phase and amplitude synchronisation in power-grid frequency fluctuations in the Nordic Grid", early access in [IEEE Access](https://doi.org/10.1109/ACCESS.2022.3150338), with doi: 10.1109/ACCESS.2022.3150338.

## Python packages
In the following follow various python scripts to obtain the figures in the scientific publication. They requires a set of python packages, which are listed below, and can easily be installed using `pip`.

Using a conventional `python >v3.4` installation, e.g., with `anaconda`, most of the standard packages should be included. These are

```code
 - numpy
 - scipy
 - pandas
 - matplotlib
```

three additional packages are needed. One for actual calculations called `MFDFA` ([MFDFA](https://github.com/LRydin/MFDFA)), which can be installed via


```code
pip install MFDFA
```

secondly, a package for plotting maps. Note: these two packages are very large and are only needed for figure no.1. If you do not want to clutter you system, we suggest you skip this (and thus disregard figure no.1). We are including here a script that plots figure no.1 without the maps, and thus does not require the heavier packages listed below

```code
pip install geopandas
pip install geoplot
```
(this could require you to install other packages via `apt-get`/`brew` and we have not tested this on Windows OS.)

## Data

Data: [Nordic_data.pkl](https://github.com/LRydin/Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid/blob/main/Nordic_data.pkl) [69.5 Mb]

The data is included here for convenience as it is used in the plots. The data is pickled and zipped, and thus in python the best way to unpickle it is with `pandas`:

```python
import pandas as pd
df = pd.read_pickle('Nordic_data.pkl', compression='zip')
```
(the `compression='zip'` is important, else pandas does not know how the data was compressed)

## Plots

 - Figure 1: [plot_1.py](https://github.com/LRydin/Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid/blob/main/plot_1.py)
 - Figure 1 (without map): [plot_1_no_map.py](https://github.com/LRydin/Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid/blob/main/plot_1_no_map.py)
 - Figure 2 (no script, made manually, included as is): [fig_2.pdf](https://github.com/LRydin/Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid/blob/main/fig_2.pdf)
 - Figure 3: [plot_3.py](https://github.com/LRydin/Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid/blob/main/plot_3.py) [the script takes a quite a bit to run]
 - Figure 4: [plot_4.py](https://github.com/LRydin/Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid/blob/main/plot_4.py) [the script as various extra plots for convenience -- it has a commented out section where the heavy calculations are done, that can take ~4 hours, thus a file with the pre-run results is also included called [pearson.npz](https://github.com/LRydin/Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid/blob/main/pearson.npz)]
 - Figure 5: [plot_4.py](https://github.com/LRydin/Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid/blob/main/plot_5.py) [the script as has a commented out section where the heavy calculations are done, that can take ~4 hours, thus a file with the pre-run results is also included called [pearson.npz](https://github.com/LRydin/Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid/blob/main/pearson.npz)]
 - Figure 6: [plot_6.py](https://github.com/LRydin/Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid/blob/main/plot_6.py) [a file with the pre-run results from the MFDFA procedure is included, called [MFDFA_data.npz](https://github.com/LRydin/Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid/blob/main/MFDFA_data.npz)]
