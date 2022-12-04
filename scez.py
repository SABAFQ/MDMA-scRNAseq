import os, sys
import pickle
import numpy as np 
import pandas as pd
from glob import glob

import scanpy as sc
import scanpy.external as sce
import anndata as ad
import scvi

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import pyplot as plot
from matplotlib.backends.backend_pdf import PdfPages

sc.settings.verbosity = 1             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100,dpi_save=300, frameon=False, figsize=(5, 5), facecolor='white')
sc.logging.print_header()

# https://stackoverflow.com/questions/21884271/warning-about-too-many-open-figures
plt.rcParams.update({'figure.max_open_warning': 0})
plt.close('all')

# https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot

SMALL_SIZE = 6
MEDIUM_SIZE = 8
BIGGER_SIZE = 10

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def add_marker_feature(adata,marker,marker_name,clusters_name,thr = 0,figsize=(10,4)):

    adata.obs[marker_name] = ''
    adata.obs.loc[adata.to_df().loc[:,marker] <= thr, marker_name] = f'{marker}-'
    adata.obs.loc[adata.to_df().loc[:,marker] > thr, marker_name] = f'{marker}+'
    
    df = pd.concat([
        adata.obs.groupby([marker_name,clusters_name]).size()[f'{marker}+'],
        adata.obs.groupby([marker_name,clusters_name]).size()[f'{marker}-']
    ],axis=1).rename(columns={0:f'{marker}+',1:f'{marker}-'})
    
    # Make some labels.
    labels = df[f'{marker}+'] / df.sum(axis=1) * 100
    labels = labels.round(decimals=1)
    labels.sort_values(ascending=False,inplace=True)
    df = df.loc[labels.index,]
    
    ax = df.plot.bar(stacked=True,rot=0,figsize=figsize)

    rects = ax.patches

    for rect, label in zip(rects, labels):
        height = rect.get_height()
        ax.text(
            rect.get_x() + rect.get_width() / 2, height + 5, str(label) + "%", 
            ha="center", va="bottom", fontsize=8
        )

    ax.set_yscale('log')
    ax.set_ylabel('# of cells')
    return ax
