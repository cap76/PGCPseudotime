import matplotlib
matplotlib.use("TkAgg")
#matplotlib.use("Agg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.path import Path
from pandas import DataFrame as df
import wishbone
import os
import sys
import platform
import pandas as pd
import tkinter as tk
import numpy as np
from tkinter import filedialog
import pickle
import random
import pylab

np.random.seed(1) 

scdata = wishbone.wb.SCData.from_csv(os.path.expanduser('~/Desktop/Code/wishbone/data/Count_data_500.csv'), 
                data_type='sc-seq', normalize=True)

#scdata = wishbone.wb.SCData.from_csv(os.path.expanduser('~/Desktop/Code/wishbone/data/sample_scseq_data_copy.csv'),
#                data_type='sc-seq', normalize=False)
scdata.run_pca()
NO_CMPNTS = 4
scdata.run_tsne(n_components=NO_CMPNTS, perplexity=20)

fig, ax = scdata.plot_gene_expression(genes = ['SOX17','DPPA3','WT1'])
#pylab.show()

scdata.run_diffusion_map()

fig, ax = scdata.plot_diffusion_components()

scdata.run_diffusion_map_correlations()
scdata.data.columns = scdata.data.columns.str.upper()
wb = wishbone.wb.Wishbone(scdata)

while True:
    try:
        wb.run_wishbone(start_cell='Cell_10', components_list=[1, 2, 3], num_waypoints=50)
    except:
         continue
    else:
         #the rest of the code
         break


fig, ax = wb.plot_wishbone_on_tsne()
vals, fig, ax = wb.plot_marker_trajectory(['WT1', 'SOX2', 'SOX17']);



pseudotime = pd.DataFrame(wb.trajectory)
branch = pd.DataFrame(wb.branch)
frames = [pseudotime,branch]

df.to_csv(pseudotime,'~/Desktop/Wishbone_pgs_soma_pseudotimes.csv')
df.to_csv(branch,'~/Desktop/Wishbone_pgs_soma_branches.csv')
