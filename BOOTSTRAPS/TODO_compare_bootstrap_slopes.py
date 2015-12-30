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
#
#
#################################
def plot_the_hist(dat,ref_point,pval,fname):
	#######################
	plt.clf()
	# fig = plt.gcf()
	hsize = 4.0
	coeff = 0.75
	fig = plt.figure(figsize=(hsize,coeff*hsize))
	heights,bins,patches = plt.hist(dat,bins=50,label='bootstrap corrs.',color='blue',linewidth=0,alpha=0.8)
	# hist height scale
	h = max(heights)
	# bin's width
	w = bins[1] - bins[0]
	#
	plt.bar(ref_point,0.2*h,w,color='red',label='observed corr.',linewidth=0)
	# plt.title("Archaeal vs Bacterial slopes correlation in bootstrap test: random Tr.Op. assignment")
	plt.xlabel("correlation")
	#
	plt.legend(loc='best',frameon=False,title="P-value:%.3f"%pval)
	plt.tight_layout()
	#
	plt.savefig(fname)
#################################
def get_pval(dat,ref_point):
	# calculate pval:
	values_above_ref = sum( cor>=ref_point for cor in dat )
	total_values = len(dat)
	return values_above_ref*1.0/total_values
#################################
#
#
#
in1 = sys.argv[1]
in2 = sys.argv[2]
ref_r = float(sys.argv[3])
out = sys.argv[4]
# # we also need a some reference value to compare with ...
# ref_r = 0.554
#################
dat1 = pd.read_csv(in1,index_col=0)
dat2 = pd.read_csv(in2,index_col=0)
#################
iters1 = dat1.shape[1]
iters2 = dat2.shape[1]
#################
correlations = []
for i in range(iters1):
    for j in range(iters2):
        x, y = dat1.icol(i), dat2.icol(j)
        a,b,r,pval,_ = stats.linregress(x,y)
        correlations.append(r)
#################
pval = get_pval(correlations,ref_r)
# plotting ...
plot_the_hist(correlations,ref_r,pval,out)
#################
#
#
#
#





