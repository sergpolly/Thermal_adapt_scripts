import os
import re
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import random as rnd
import matplotlib as mpl

aacids  = list("CMFILVWYAGTSNQDEHRKP")
freq_keys = sorted(aacids)


# mpl.rcParams['font.family'] = 'sans-serif'
font = {'family' : 'sans-serif',
        #'weight' : 'bold',
        'size'   : 8}

mpl.rc('font', **font)

# slopes to compare with  ...
fnames = {}

in1 = sys.argv[1]
in2 = sys.argv[2]


dat1 = pd.read_csv(in1,index_col=0)
dat2 = pd.read_csv(in2,index_col=0)

iters1 = dat1.shape[1]
iters2 = dat2.shape[1]


correlations = []
for i in range(iters1):
    for j in range(iters2):
        x, y = dat1.icol(i), dat2.icol(j)
        a,b,r,pval,_ = stats.linregress(x,y)
        correlations.append(r)



