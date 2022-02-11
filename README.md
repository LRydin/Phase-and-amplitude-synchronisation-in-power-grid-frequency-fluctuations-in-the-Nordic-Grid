# Phase-and-amplitude-synchronisation-in-power-grid-frequency-fluctuations-in-the-Nordic-Grid
Related python code to the scientific publication "Phase and amplitude synchronisation in power-grid frequency fluctuations in the Nordic Grid", early access in [IEEE Access](https://doi.org/10.1109/ACCESS.2022.3150338), with doi: 10.1109/ACCESS.2022.3150338.

In the following follow various python scripts to obtain the figures in the scientific publication. They requires a set of python packages, which are listed below, and can easily be installed using `pip`.

Using a conventional `python >v3.4` installation, e.g., with `anaconda`, most of the standard packages should be included. These are

```code
 - numpy
 - scipy
 - pandas
 - matplotlib
```

three additional packages are needed. One for the actual calculations called `MFDFA` ([MFDFA](https://github.com/LRydin/MFDFA)), which can be installed via


```code
pip install MFDFA
```

secondly, a package for plotting maps. Note: these two packages are very large and are only needed for figure no.1. If you do not want to clutter you system, we suggest you skip this (and thus disregard figure no.1). We are including here a script that plots figure no.1 without the maps, and thus does not require the heavier packages listed below

```code
pip install geopandas
pip install geoplot
```
(this could require you to install other packages via `apt-get`/`brew` and we have not tested this on Windows OS.)
