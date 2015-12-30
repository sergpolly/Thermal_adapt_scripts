import re
import os
import sys
import pandas as pd
import numpy as np
from scipy import stats as st

# for plotting ...
import matplotlib.pyplot as plt
import matplotlib as mpl
#
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
font = {#'family' : 'sans-serif',
        #'weight' : 'bold',
        'size'   :9}
rc('font', **font)
# # data loading ...

aacids = list('ACDEFGHIKLMNPQRSTVWY')

def plot_comparison(ax,x,y,xlab='x',ylab='y'):
    #
    ax.plot(x,y,'bo',visible=False)
    for i,aa in enumerate(aacids):
        ax.text(x[i], y[i], r"\texttt{%s}"%aa, transform=ax.transData,color='dimgray',horizontalalignment='center',verticalalignment='center',fontsize=12,fontweight='bold')
    ###################################################
    a,b,r,pval,_ = st.linregress(x,y)
    #
    middle = 0.5*(x.min()+x.max())
    delta = 1.1*(x.max()-x.min())
    x_range = [middle-0.5*delta,middle+0.5*delta]
    #############################################
    fff = np.vectorize(lambda x: a*x + b)
    lin_fit, = ax.plot(x_range,fff(x_range),color='blue',linewidth=1.5,linestyle='-',zorder=100,label='linear fit, R=%.3f, p=%.3f'%(r,pval))
    ###################################################
    ###################################################
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ###################################################
    ###################################################
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.legend(loc='best',frameon=False)
    # ax.legend((lin_fit,),('linear fit, R=%.3f, p=%.3f'%(r,pval)),loc='best',frameon=False)
    plt.tight_layout(pad=0.4, h_pad=None, w_pad=None)
    ###################################################
    ###################################################
    # axis limits ...
    middle = 0.5*(x.min()+x.max())
    delta = 1.3*(x.max()-x.min())
    x_lims = [middle-0.5*delta,middle+0.5*delta]
    middle = 0.5*(y.min()+y.max())
    delta = 1.3*(y.max()-y.min())
    y_lims = [middle-0.5*delta,middle+0.5*delta]
    # ###################################################
    ax.set_xlim(x_lims)
    ax.set_ylim(y_lims)
    ###################################################
    # ax.locator_params(axis='both',tight=True)#nbins=10)
    # ax.locator_params(axis='x',nbins=5)
    ###################################################



def get_axis(xbins = 2, ybins = 3,vcoeff=1.0):
    plt.clf()
    # create grid of subplots on a large figure canvas
    # share x&y axes among all plots
    fig, ax = plt.subplots(ybins, xbins, figsize=(7.5,ybins*vcoeff*7.5/xbins))#, sharex=True, sharey=True)
    # # no space between subplots
    # fig.subplots_adjust(hspace=0.0, wspace=0.08)
    # # make some room for the axes' labels
    # l,b,r,t = 0.05,0.12,0.98,0.92
    # w, h = r-l, t-b
    # fig.subplots_adjust(bottom=b, left=l, right=r, top=t)
    # #
    # # returning axis ...
    return (fig,ax)


def get_arch_bact_slopes(combination):
    # getting file name based on criteria combination ...
    def get_exp_fname(cds_criteria,org_criteria,kingdom):
        if kingdom in ['arch','archaea']:
            return "exp_MTAD_%s_%s.arch.summary"%(cds_criteria,org_criteria)
        elif kingdom in ['bact','bacteria']:
            return "exp_MTAD_%s_%s.bact.summary"%(cds_criteria,org_criteria)
        else:
            raise TypeError('only archaeal and bacterial kingdoms are supported!')
    #
    #
    # get archaeal slopes ...
    arch_fname = get_exp_fname(combination['cds_criteria'],
                        combination['org_criteria'],
                        'arch')
    arch = pd.read_csv(arch_fname,index_col=0)
    # get bacterial slopes ...
    bact_fname = get_exp_fname(combination['cds_criteria'],
                        combination['org_criteria'],
                        'bact')
    bact = pd.read_csv(bact_fname,index_col=0)
    #
    return (arch,bact)


def get_axis_labels(combination):
    # unpack criteria ...
    cds_crit,org_crit = combination['cds_criteria'], combination['org_criteria']
    label = ""
    if cds_crit == 'cai':
        label += "top 10\% CAI"
    elif cds_crit == 'ribo':
        label += "ribosomal proteins"
    elif cds_crit == 'cai_noribo':
        label += "top 10\% CAI excl. ribosomal"
    elif cds_crit == 'all':
        label += "proteome"
    else:
        raise ValueError("cds criteria not supported!")
    #
    if org_crit == 'all':
        arch_label = ' '.join([label,r"(archaea), 1/\textdegree C"])
        bact_label = ' '.join([label,r"(bacteria), 1/\textdegree C"])
    elif org_crit == 'trop':
        arch_label = ' '.join([label,r"(CUS archaea), 1/\textdegree C"])
        bact_label = ' '.join([label,r"(CUS bacteria), 1/\textdegree C"])
    #
    return (arch_label,bact_label)


cds_crits = ['cai','ribo','cai_noribo','all']
org_crits = ['trop','all']

####################
# THERE ARE JUST 4 COMBINATIONS OF CRITS, SO JUST DO THEM MANUALLY ...
combinations = [{'cds_criteria':'all','org_criteria':'all'},
                {'cds_criteria':'ribo','org_criteria':'all'},
                {'cds_criteria':'all','org_criteria':'trop'},
                {'cds_criteria':'cai','org_criteria':'trop'}]
###########
#
# combinations = [(cc,oo) for cc in cds_crits for oo in org_crits]
# combinations = [[{'cds_criteria':cc,'org_criteria':oo} for cc in cds_crits] for oo in org_crits]
# # # with trop
# # combinations[0][0,2,3]
# # # no trop ...
# # combinations[1][0,1,3]


xsize, ysize = 2, 2
# # draw the thing ...
fig,ax = get_axis(xsize, ysize, vcoeff = 0.8)

for cid,criteria_combination in enumerate(combinations):
    # axis coordinates ...
    i, j = cid/xsize, cid%xsize
    #
    arch,bact = get_arch_bact_slopes(criteria_combination)
    alab,blab = get_axis_labels(criteria_combination)
    plot_comparison(ax[i,j],arch['exp_D'],bact['exp_D'],xlab=alab,ylab=blab)





# for i,(for_trop,for_notrop) in enumerate(zip([0,2,3],[0,1,3])):
#     #######################
#     criteria_combination = combinations[0][for_trop]
#     print i,criteria_combination
#     arch,bact = get_arch_bact_slopes(criteria_combination)
#     alab,blab = get_axis_labels(criteria_combination)
#     plot_comparison(ax[i,0],arch['exp_D'],bact['exp_D'],xlab=alab,ylab=blab)
#     #######################
#     criteria_combination = combinations[1][for_notrop]
#     print i,criteria_combination
#     arch,bact = get_arch_bact_slopes(criteria_combination)
#     alab,blab = get_axis_labels(criteria_combination)
#     plot_comparison(ax[i,1],arch['exp_D'],bact['exp_D'],xlab=alab,ylab=blab)



plt.savefig("figure2.pdf")

































