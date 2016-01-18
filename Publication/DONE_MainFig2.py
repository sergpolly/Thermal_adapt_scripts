import re
import os
import sys
import pandas as pd
import numpy as np
from scipy import stats as st


# for plotting ...
import subprocess as sub
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

def plot_comparison(ax,x,y,xlab='x',ylab='y',xlims=None,ylims=None):
    #
    ax.plot(x,y,'bo',visible=False)
    for i,aa in enumerate(aacids):
        ax.text(x[i], y[i], r"\texttt{%s}"%aa, transform=ax.transData,color='black',horizontalalignment='center',verticalalignment='center',fontsize=13,fontweight='bold')
    ###################################################
    a,b,r,pval,_ = st.linregress(x,y)
    #
    middle = 0.5*(x.min()+x.max())
    delta = 1.1*(x.max()-x.min())
    x_range = [middle-0.5*delta,middle+0.5*delta]
    #############################################
    label_f = lambda r,pval:("linear fit, R=%.3f, "%r) + (r"$p=%.3f$"%pval if pval>=0.001 else r"$p<0.001$")
    fff = np.vectorize(lambda x: a*x + b)
    lin_fit, = ax.plot(x_range,fff(x_range),color='silver',linewidth=1.5,linestyle='-',zorder=100,label=label_f(r,pval))
    y_eq_x, = ax.plot(x_range,x_range,color='red',linewidth=1.0,linestyle='--',zorder=102, label=r"y=x guide line")
    ###################################################
    ###################################################
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ###################################################
    ###################################################
    ax.set_xlabel(xlab)#,labelpad=1.5)
    ax.set_ylabel(ylab)#,labelpad=0.5)
    ax.legend(loc='upper left',frameon=False,handlelength=1.5,handletextpad=0.1)
    # ax.legend((lin_fit,),('linear fit, R=%.3f, p=%.3f'%(r,pval)),loc='best',frameon=False)
    plt.tight_layout(pad=0.4, h_pad=None, w_pad=None)
    ###################################################
    ###################################################
    if (xlims is None) or (ylims is None):
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
    else:
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
    ###################################################
    # ax.locator_params(axis='both',tight=True)#nbins=10)
    # ax.locator_params(axis='x',nbins=5)
    ###################################################



def get_axis(xbins = 2, ybins = 3,vcoeff=1.0):
    plt.clf()
    # create grid of subplots on a large figure canvas
    # share x&y axes among all plots
    # fig, ax = plt.subplots(ybins, xbins, figsize=(7.5,ybins*vcoeff*7.5/xbins))#, sharex=True, sharey=True)
    x_fig_size = 3.34
    fig, ax = plt.subplots(ybins, xbins, figsize=(x_fig_size,vcoeff*x_fig_size))#, sharex=True, sharey=True)
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
    #
    label = ""
    if cds_crit == 'cai':
        label += " (top 10\% CAI)"
    elif cds_crit == 'ribo':
        label += " (ribosomal proteins)"
    elif cds_crit == 'cai_noribo':
        label += " (top 10\% CAI excl. ribosomal)"
    elif cds_crit == 'all':
        label += " (proteome)"
    else:
        raise ValueError("cds criteria not supported!")
    #
    if org_crit == 'all':
        arch_label = ''.join(["archaea slopes",label,r", 1/\textdegree C"])
        bact_label = ''.join(["bacteria slopes",label,r", 1/\textdegree C"])
    elif org_crit == 'trop':
        arch_label = ''.join(["CUS archaea slopes",label,r", 1/\textdegree C"])
        bact_label = ''.join(["CUS bacteria slopes",label,r", 1/\textdegree C"])
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


xsize, ysize = 1,1
xlims = [(-0.035,0.035),(-0.035,0.035),(-0.035,0.035),(-0.035,0.035)]
ylims = [(-0.065,0.075),(-0.055,0.055),(-0.055,0.065),(-0.04,0.05)]
# # draw the thing ...
for cid,criteria_combination in enumerate(combinations):
    fig,ax = get_axis(xsize, ysize, vcoeff = 0.85)
    # # axis coordinates ...
    # i, j = cid/xsize, cid%xsize
    #
    arch,bact = get_arch_bact_slopes(criteria_combination)
    alab,blab = get_axis_labels(criteria_combination)
    plot_comparison(ax,arch['exp_D'],bact['exp_D'],xlab=alab,ylab=blab,xlims=xlims[cid],ylims=ylims[cid])
    plt.savefig("figure2_%s.png"%'abcd'[cid],dpi=300)

cmd = "montage "+\
"figure2_a.png "+\
"figure2_b.png "+\
"figure2_c.png "+\
"figure2_d.png "+\
"-geometry +0+0 "+\
"-tile 2x "+\
"-border 0 "+\
"-bordercolor white "+\
"-background white "+\
"figure2.png"
# "-geometry +0-15 "+\

retcode = sub.call(cmd,shell=True)








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




































