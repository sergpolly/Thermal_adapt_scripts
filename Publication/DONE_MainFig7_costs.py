import os
import sys
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats as st
import random as rnd
import numpy as np
#
#
#
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullFormatter
import scipy.interpolate as interpol
# font = {'family' : 'sans-serif',
#         #'weight' : 'bold',
#         'size'   :9}
# #
# mpl.rc('font', **font)
# #
# #
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# #
#
mpl.rcParams['text.latex.preamble'] = [
       r'\usepackage{textcomp}',   # i need upright \micro symbols, but you need...
       # r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]
#
#
font = {#'family' : 'sans-serif',
        #'weight' : 'bold',
        'size'   :9}
rc('font', **font)
# # data loading ...
#
#
#
#
aacids = list('CMFILVWYAGTSNQDEHRKP')
aa_combinations = ['IVYWREL', 'DEKR', 'AGNQSTHY', 'MPCLVWIF', 'ILVM']
#
#
root_path = os.path.expanduser('~')
bact_path = os.path.join(root_path,'GENOMES_BACTER_RELEASE69/genbank')
arch_path = os.path.join(root_path,'GENOMES_ARCH_SEP2015')
# SOME ARCHAEAL DATA ...
arch        = pd.read_csv(os.path.join(arch_path,'summary_organisms_interest.dat'))
arch_nohalo = pd.read_csv(os.path.join(arch_path,'summary_organisms_interest_no_halop.dat'))


# SOME BACTERIAL DATA ...
# complete genomes only ...
bact        = pd.read_csv(os.path.join(bact_path,'env_catalog_compgenome.dat'))

#
# bacter proteomic summary ...
bact_prot = pd.read_csv(os.path.join(bact_path,'proteome_all.dat'))
#
# arch proteomic summary ...
arch_prot = pd.read_csv(os.path.join(arch_path,'proteome_arch.dat'))
arch_prot[aacids] = arch_prot[aacids]*100.0

arch_dat = pd.merge(arch,arch_prot,on='assembly_accession')
arch_nohalo_dat = pd.merge(arch_nohalo,arch_prot,on='assembly_accession')
arch_halo_dat = arch_dat[~arch_dat['assembly_accession'].isin(arch_nohalo['assembly_accession'])]

bact_dat = pd.merge(bact,bact_prot,on='GenomicID')
bact_dat[aacids] = bact_dat[aacids]*100.0


##################################
###    COST VECTORS LOADING    ###
##################################
# cost vectors loading ...
cost_vec_path = '.'
akashi = os.path.join(cost_vec_path,'akashi-cost.d')
argentina = os.path.join(cost_vec_path,'argentina-cost.d')
##################
akashi_cost = pd.read_csv(akashi,header=None,sep=' ')
argentina_cost = pd.read_csv(argentina,header=None,sep=' ')
##################
akashi_cost.set_index(0,inplace=True)
argentina_cost.set_index(0,inplace=True)
# loaded ...
#######################################################################################################


##########################################
# what data to substitute to the screipt ...
if len(sys.argv)<=1:
    raise ValueError("arch ot bact? choose the kingdom to proceed!")
kingdom = sys.argv[1]
if kingdom == 'bact':
    dat = bact_dat
elif kingdom == 'arch':
    dat = arch_nohalo_dat



def label(rr,pp):
    if (pp<0.0008):
        label = '$R=%.2f,\ p<0.001$'%rr
    else:
        label = '$R=%.2f,\ p=%.3f$'%(rr,pp)
    return label



def get_lims(dat,coeff=1.0):
    lims = (dat.min(),dat.max())
    middle = 0.5*(lims[0] + lims[1])
    half = abs(lims[1] - middle) # = abs(lims[0] - middle)
    return (middle-coeff*half,middle+coeff*half)


topt = 'OptimumTemperature'
vmin,vmax = get_lims(dat['GC'],coeff=1.0)
xmin,xmax = get_lims(dat[topt],coeff=1.03)
regress_lims = get_lims(dat[topt],coeff=0.95)

####################################################
# plot the simulated proteome's cost (both Akashi and Argentina)
plt.clf()
x_fig_size = 7.4
v_coeff = 0.4
fig = plt.figure(figsize=(x_fig_size,v_coeff*x_fig_size))
# between axes ...
hor_space = 0.07
# axes info ...
left = 0.09
bottom = 0.15
width = 0.5*(0.87 - left - hor_space)
height = 0.975 - bottom
# bottom axes 
ax_left = plt.axes([left, bottom, width, height])
left += (width + hor_space)
# top axes 
ax_right = plt.axes([left, bottom, width, height])
#
#
print "Akashi corrs (a,b,R,P) ..."
#
evolved_proteome_akashi_cost = (dat[aacids]/100.0).dot(akashi_cost)
evolved_proteome_akashi_cost = evolved_proteome_akashi_cost[evolved_proteome_akashi_cost.columns[0]]
#
ymin,ymax = get_lims(evolved_proteome_akashi_cost,coeff=1.1)
#
scatter = ax_left.scatter(dat[topt],evolved_proteome_akashi_cost,s=65,c=dat['GC'],edgecolor='none',vmin=vmin,vmax=vmax,cmap=plt.get_cmap('jet'))
# linear regression & fit ...
a,b,r,pval,_ = st.linregress(dat[topt],evolved_proteome_akashi_cost)
print "all:    ",a,b,r,pval
t_range = np.asarray(regress_lims)
ax_left.plot(t_range,a*t_range+b,'-',color='dimgray',lw=2,label=label(r,pval))
ax_left.legend(loc='best',frameon=False)
#
ax_left.set_xlabel(r'OGT,\textdegree C')
ax_left.set_ylabel('AA synthesis cost, ATP')
#
#################################
ax_left.yaxis.set_ticks_position('left')
ax_left.xaxis.set_ticks_position('bottom')
#
ax_left.set_xlim((xmin,xmax))
ax_left.set_ylim((ymin,ymax))
#################################
print "Argentina corrs (a,b,R,P) ..."
#
evolved_proteome_argentina_cost = (dat[aacids]/100.0).dot(argentina_cost)
evolved_proteome_argentina_cost = evolved_proteome_argentina_cost[evolved_proteome_argentina_cost.columns[0]]
#
ymin,ymax = get_lims(evolved_proteome_argentina_cost,coeff=1.1)
#
scatter = ax_right.scatter(dat[topt],evolved_proteome_argentina_cost,s=65,c=dat['GC'],edgecolor='none',vmin=vmin,vmax=vmax,cmap=plt.get_cmap('jet'))
# linear regression & fit ...
a,b,r,pval,_ = st.linregress(dat[topt],evolved_proteome_argentina_cost)
print "all:    ",a,b,r,pval
t_range = np.asarray(regress_lims)
ax_right.plot(t_range,a*t_range+b,'-',color='dimgray',lw=2,label=label(r,pval))
ax_right.legend(loc='best',frameon=False)
#
ax_right.set_xlabel(r'OGT,\textdegree C')
ax_right.set_ylabel('AA maintenance cost, ATP/time')
#################################
ax_right.yaxis.set_ticks_position('left')
ax_right.xaxis.set_ticks_position('bottom')
#
ax_right.set_xlim((xmin,xmax))
ax_right.set_ylim((ymin,ymax))
#
cax = fig.add_axes([left+width+0.01,bottom,0.03,height])
cbar = fig.colorbar(scatter,cax=cax,orientation='vertical')
ticks = 5*(pd.np.arange(vmin//5,vmax//5)+1)
ticklabels = map(str,ticks)
cbar.set_ticks(ticks)
cbar.set_ticklabels(ticklabels)
cbar.set_label('GC content, \%')
#
# fig.savefig(os.path.join(results_path,"%s.png"%fname),dpi=300)
fig.savefig("Fig7.%s.pdf"%kingdom)









































