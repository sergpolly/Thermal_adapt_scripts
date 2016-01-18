
import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess as sub
import pandas as pd
import numpy as np
import sys
import re
import os
import scipy.stats as st

from matplotlib.ticker import MaxNLocator
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Rectangle

import scipy.interpolate as interpol

import matplotlib.patches as patches
import matplotlib.transforms as transforms
#
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
#
font = {#'family' : 'sans-serif',
        #'weight' : 'bold',
        'size'   :9}
rc('font', **font)
# # data loading ...



# font = {'family' : 'sans-serif',
#         # 'sans-serif':['Helvetica'],
#         #'weight' : 'bold',
#         'size'   : 7.5}

# mpl.rc('font', **font)

# constants ...
# ALPHABET=20
# amino acis alphabet ...
# aacids = data.columns.values[:ALPHABET]

aacids = sorted('CMFILVWYAGTSNQDEHRKP')
# colors ...
# cols = [plt.cm.gist_rainbow(ccc) for ccc in  np.linspace(0.0,1.0,num=20)]

#############
#############################################################################
kelly_colors_hex = [
    u'#FFB300', # Vivid Yellow
    u'#803E75', # Strong Purple
    u'#FF6800', # Vivid Orange
    u'#A6BDD7', # Very Light Blue
    u'#C10020', # Vivid Red
    u'#CEA262', # Grayish Yellow
    u'#817066', # Medium Gray

    u'#007D34', # Vivid Green
    u'#F6768E', # Strong Purplish Pink
    u'#00538A', # Strong Blue
    u'#FF7A5C', # Strong Yellowish Pink
    u'#53377A', # Strong Violet
    u'#FF8E00', # Vivid Orange Yellow
    u'#B32851', # Strong Purplish Red
    u'#F4C800', # Vivid Greenish Yellow
    u'#7F180D', # Strong Reddish Brown
    u'#93AA00', # Vivid Yellowish Green
    u'#593315', # Deep Yellowish Brown
    u'#F13A13', # Vivid Reddish Orange
    u'#232C16' # Dark Olive Green
    ]
#############################################################################
cols = kelly_colors_hex
#############



def get_axis(xbins = 2, ybins = 3, yext=1.0):
    plt.clf()
    # create grid of subplots on a large figure canvas
    # share x&y axes among all plots
    fig, ax = plt.subplots(ybins, xbins, figsize=(7.5,ybins*yext*7.5/xbins))#, sharex=True, sharey=True)
    # # no space between subplots
    fig.subplots_adjust(hspace=0.05, wspace=0.2)
    # # make some room for the axes' labels
    l,b,r,t = 0.08,0.08,0.98,0.97
    # w, h = r-l, t-b
    fig.subplots_adjust(bottom=b, left=l, right=r, top=t)
    # #
    # # returning axis ...
    return (fig,ax)


# getting file name based on criteria combination ...
def get_exp_fname(cds_criteria,org_criteria,kingdom):
    if kingdom in ['arch','archaea']:
        return "exp_MTAD_%s_%s.arch.summary"%(cds_criteria,org_criteria)
    elif kingdom in ['bact','bacteria']:
        return "exp_MTAD_%s_%s.bact.summary"%(cds_criteria,org_criteria)
    else:
        raise TypeError('only archaeal and bacterial kingdoms are supported!')


def get_slopes_comparison_labels(combination,kingdom):
    # unpack criteria ...
    cds_crit,org_crit = combination['cds_criteria'], combination['org_criteria']
    label = r"observed slopes "
    ######################
    def translate_kingdom(kingdom):
        if kingdom in ['arch','archaea']:
            return 'archaea'
        elif kingdom in ['bact','bacteria']:
            return 'bacteria'
        else:
            raise ValueError("Only archaeal and bacterial kingdoms are supported!")
    ######################
    the_kingdom = translate_kingdom(kingdom)
    ######################
    if org_crit == 'all':
        if cds_crit == 'cai':
            label += r"(top 10\% CAI %s), 1/\textdegree C"%the_kingdom
        elif cds_crit == 'ribo':
            label += r"(ribosomal proteins %s), 1/\textdegree C"%the_kingdom
        elif cds_crit == 'cai_noribo':
            label += r"(top 10\% CAI excl. ribosomal %s), 1/\textdegree C"%the_kingdom
        elif cds_crit == 'all':
            label += r"(proteome %s), 1/\textdegree C"%the_kingdom
        else:
            raise ValueError("cds criteria not supported!")
    elif org_crit == 'trop':
        if cds_crit == 'cai':
            label += r"(top 10\%" + r" CAI, CUS %s), 1/\textdegree C"%the_kingdom
        elif cds_crit == 'ribo':
            label += r"(ribosomal proteins, CUS %s), 1/\textdegree C"%the_kingdom
        elif cds_crit == 'cai_noribo':
            label += r"(top 10\%" + r" CAI excl. ribosomal, CUS %s), 1/\textdegree C"%the_kingdom
        elif cds_crit == 'all':
            label += r"(proteome, CUS %s), 1/\textdegree C"%the_kingdom
        else:
            raise ValueError("cds criteria not supported!")
    #######################################
    return label


# cds_crits = ['cai','ribo','cai_noribo','all']
# org_crits = ['trop','all']
# combinations = [(cc,oo) for cc in cds_crits for oo in org_crits]
# combinations = [[{'cds_criteria':cc,'org_criteria':oo} for cc in cds_crits] for oo in org_crits]



#######################################################################################################
# TODO  GENERATE APPROPROATE FILES HERE ...
if len(sys.argv)<=3:
    raise TypeError("Use command line argument to enter CDS and organismal criteria!")
else:
    the_combination = {'cds_criteria':sys.argv[1],'org_criteria':sys.argv[2]}
    the_kingdom = sys.argv[3]
#
exp_fname = get_exp_fname(the_combination['cds_criteria'],the_combination['org_criteria'],the_kingdom)
#
# "observed slopes (proteome bacteria), $1/^{o}C$"
slopes_label = get_slopes_comparison_labels(the_combination,the_kingdom)
#
#
#
#
data_exp = pd.read_csv(exp_fname,index_col=0)
# loaded ...
# give em names, as I figured, indexes and column names are awesome things ...
exp_T = data_exp['exp_T']
exp_M = data_exp['exp_M']
exp_A = data_exp['exp_A']
exp_D = data_exp['exp_D']
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

################################################
# WHERE RESULTS(PICTURES) SHOUDL GO ...
################################################
results_path = '.'
################################################

#####################################################
# SIMULATION DATA ...
# PATH to simulations data ...
# protein_design/the_simulations
simul_path = os.path.join(os.path.expanduser('~'),
                            "Dropbox (UMASS MED - BIB)",
                            "protein_design",
                            "the_simulations",
                            "Correct_MJ99_Argentina_PUBLICATION")
# data files and their names ...
simul_fname = os.path.join(simul_path,"simulation.dat")
# loading simul data ...
data_sorted = pd.read_csv(simul_fname)
#####################################################


####################################################
# SOME PLOTTING FUNCTIONS ...
####################################################
def plot_broom(MFFTW_case, data):
    # plt.clf()
    fig = plt.figure(figsize=(5.5,4.5))
    # ax = fig.add_subplot(1,1,1)
    ax = fig.add_axes((0.09,0.09,0.78,0.85))
    #
    # rect = l,b,w,h fig.add_axes(rect)
    #
    M,FFT,W=MFFTW_case
    sub_data = data[(data.Model == M)&(data.W_coeff == W)&(data.FFTemp == FFT)]
    for aa, col in zip(aacids,cols):
        ax.plot(sub_data.Temp,sub_data[aa]*100.0,'o-',ms=4, lw=2, mew=0.0, color=col, label=aa)
    #
    # leg = plt.legend(loc='best',fontsize=13)
    leg = ax.legend(loc='center left', bbox_to_anchor=(0.999, 0.5),fontsize=9)
    # leg = plt.legend(loc='best')
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    # wider markers in the legend ...
    for legobj in leg.legendHandles:
        legobj.set_linewidth(6.0)
    #
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    #
    ax.set_xlabel('Temperature, units')
    ax.set_ylabel('composition, %')
    #
    # ax.set_xlim((0.0,6.1))
    # ax.set_ylim((2.5,10.0))
    #
    # fig.tight_layout()
    #
    plt.title(" case %d %.4f FFT=%s" % (M,W,FFT))
    plt.savefig(os.path.join(results_path,"der_case_%d_%.4f_FFT%s.ghpcc.png"%(M,W,FFT)),dpi=300)
###############################
##########################
##########################
##########################
###############################
def plot_broom_groups(MFFTW_case, data):
    # plt.clf()
    fig = plt.figure(figsize=(5.5,4.5))
    # ax = fig.add_subplot(1,1,1)
    ax = fig.add_axes((0.09,0.09,0.85,0.85))
    #
    # rect = l,b,w,h fig.add_axes(rect)
    #
    M,FFT,W=MFFTW_case
    sub_data = data[(data.Model == M)&(data.W_coeff == W)&(data.FFTemp == FFT)]
    charged = sub_data[['D','E','K','R']].sum(axis=1)
    hydrophob = sub_data[['A', 'G', 'N', 'Q', 'S', 'T', 'H', 'Y']].sum(axis=1)
    hydrophil = sub_data[['M', 'P', 'C', 'L', 'V', 'W', 'I', 'F']].sum(axis=1)    
    IVYWREL = sub_data[['I', 'V', 'Y', 'W', 'R', 'E', 'L']].sum(axis=1)    
    ax.plot(sub_data.Temp,charged*100.0,'o-',ms=6, lw=4, mew=0.0, color='r', label='DEKR')
    ax.plot(sub_data.Temp,hydrophob*100.0,'o-',ms=6, lw=4, mew=0.0, color='g', label='AGNQSTHY')
    ax.plot(sub_data.Temp,hydrophil*100.0,'o-',ms=6, lw=4, mew=0.0, color='b', label='MPCLVWIF')
    ax.plot(sub_data.Temp,IVYWREL*100.0,'o-',ms=6, lw=4, mew=0.0, color='m', label='IVYWREL')
    # for aa, col in zip(aacids,cols):
    #     ax.plot(sub_data.Temp,sub_data[aa]*100.0,'o-',ms=4, lw=2, mew=0.0, color=col, label=aa)
    #
    leg = plt.legend(loc='center left',fontsize=15)
    # leg = ax.legend(loc='center left', bbox_to_anchor=(0.999, 0.5),fontsize=9)
    # leg = plt.legend(loc='best')
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_edgecolor('white')
    # wider markers in the legend ...
    for legobj in leg.legendHandles:
        legobj.set_linewidth(6.0)
    #
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    #
    ax.set_xlabel('Temperature, units')
    ax.set_ylabel('composition, %')
    #
    # ax.set_xlim((0.0,6.1))
    # ax.set_ylim((2.5,10.0))
    #
    # fig.tight_layout()
    #
    plt.title(" case %d %.4f FFT=%s"%(M,W,FFT))
    plt.savefig(os.path.join(results_path,"groups_case_%d_%.4f_FFT%s.ghpcc.png"%(M,W,FFT)),dpi=300)
##########################
##########################
##########################
##########################
##########################
def get_derivative_data_new(data,dT):
    sub_data = data.copy()
    der_data = data.copy()
    for aa in aacids:
        # new simple ...
        der_data[aa].iloc[0:-1] = (sub_data[aa].iloc[1:].values-sub_data[aa].iloc[0:-1].values)/dT
        der_data[aa].iloc[-1] = (sub_data[aa].iloc[-1]-sub_data[aa].iloc[-2])/dT
    return der_data
###########################
###########################
def get_derivative_data_super(data,dT):
    result = (data[aacids].iloc[1:].reset_index(drop=True)-data[aacids].iloc[0:-1].reset_index(drop=True))
    result = result.append(result.iloc[-1])
    result = result/float(dT)
    result = result.set_index(data.index)
    return result
###########################
####################################################
# SOME PLOTTING FUNCTIONS ...
####################################################





###########################################
###########################################
###########################################
# T & M data preparation ... correlations are being calculated ...
###########################################
###########################################
###########################################
mod = 64
fft = 'MJ99'
#
cors = []
# wopMT = 0.06
# [0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15]
for www in  data_sorted[(data_sorted.Model==mod)&(data_sorted.FFTemp==fft)].W_coeff.unique():
    data_MT = data_sorted[(data_sorted.Model==mod)&(data_sorted.FFTemp==fft)&(data_sorted.W_coeff==www)]
    data_MT.set_index('Temp',verify_integrity=True,inplace=True)
    #
    # print "getting to the place ...."
    ###############################################
    ###############################################
    ###############################################
    dT = data_MT.index[1] - data_MT.index[0]
    # print "place  1 ...."
    # print data_MT
    der_dat = get_derivative_data_super(data_MT,float(dT))
    der_dat = der_dat.append(exp_D)
    res_corr_D = der_dat[aacids].transpose().corr(method='pearson')
    no_self_corr_D = (res_corr_D.index!='exp_D')
    cors_D = res_corr_D['exp_D'].loc[no_self_corr_D]    
    ###############################################
    ###############################################
    ###############################################
    #
    # break
    #
    # add 2 exp columns to use .corr afterwards ...
    data_MT = data_MT.append(exp_T)
    data_MT = data_MT.append(exp_M)
    data_MT = data_MT.append(exp_A)
    # correlation method 'pearson' or 'spearman' ...
    results_corr = data_MT[aacids].transpose().corr(method='pearson')
    # we need corrs of expM and expT with simulated data, not with each other ...
    no_self_corr = (results_corr.index!='exp_M')&(results_corr.index!='exp_T')&(results_corr.index!='exp_A')
    cors_M = results_corr['exp_M'].loc[no_self_corr]
    cors_T = results_corr['exp_T'].loc[no_self_corr]
    cors_A = results_corr['exp_A'].loc[no_self_corr]
    # form the DataFrame with TM correlations as columns and corresponding wcf as well ...
    cors_df = pd.DataFrame({'exp_M':cors_M,'exp_T':cors_T,'exp_A':cors_A,'exp_D':cors_D,'W_coeff':[www,]*cors_M.size})
    cors_df.index.name = 'Temp'
    cors.append(cors_df)
    #########################
    # T_M = pd.np.asarray(data_MT['Temp'])[corsM.argmax()]
    # T_T = pd.np.asarray(data_MT['Temp'])[corsT.argmax()]


# concat all DataFrames corresponding to different wcf ...
cors_TM = pd.concat(cors).reset_index()

# so far, these are 2 different w and T, but w are close to each other ...
wopM = cors_TM.loc[cors_TM.exp_M.argmax()]['W_coeff']
wopT = cors_TM.loc[cors_TM.exp_T.argmax()]['W_coeff']
T_M = cors_TM.loc[cors_TM.exp_M.argmax()]['Temp']
T_T = cors_TM.loc[cors_TM.exp_T.argmax()]['Temp']


# now we should plot these things


##############################
# data for the next figure ...
RMT_max = [(www,cors_wTM['exp_M'].max(),cors_wTM['exp_T'].max()) for www,cors_wTM in cors_TM.groupby('W_coeff')]
RMT_max = pd.DataFrame(RMT_max,columns=['W_coeff','RM_max','RT_max'])


# combination ...
# RMT_max['RMT_max'] = RMT_max['RT_max']+RMT_max['RM_max']
RMT_max['RMT_max'] = RMT_max['RT_max']/RMT_max['RT_max'].max()+RMT_max['RM_max']/RMT_max['RM_max'].max()
wopMT = RMT_max['W_coeff'][RMT_max['RMT_max'].argmax()]


# print "Optimal temperatures for M,T: %.2f,%.2f reached at wop=%.3f"%(T_M,T_T,wopMT)
print "Optimal parameters for M,T:  T_M=%.2f at wM=%.2f  and  T_T=%.2f at wT=%.2f"%(T_M,wopM,T_T,wopT)
print "Overall optimal balance coefficient wMT=%.2f"%wopMT
cors_TMop = cors_TM[cors_TM.W_coeff==wopMT]
Top_M=cors_TMop.loc[cors_TMop.exp_M.argmax()].Temp
Top_T=cors_TMop.loc[cors_TMop.exp_T.argmax()].Temp
print "Overall optimal temperature range at wMT=%.2f is (%.2f,%.2f)"%(wopMT,Top_M,Top_T)
# # borderpad          the fractional whitespace inside the legend border
# # labelspacing       the vertical space between the legend entries
# # handlelength       the length of the legend handles
# # handletextpad      the pad between the legend handle and text
# # borderaxespad      the pad between the axes and legend border
# # columnspacing      the spacing between columns

# calculations are over, only plotting down below ...
###########################################################################



fig_TMAD,ax_TMAD = get_axis(xbins = 2, ybins = 2, yext =0.85)
#################################
# first figure would be T & M - thing ...
# bottom axes 
ax_bottom = ax_TMAD[1,0]
# top axes 
ax_top = ax_TMAD[0,0]
################################
################################
################################
wcf_plt_M = [0.,0.01,0.03]+[wopM,]+[0.11,0.15]
wcf_plt_T = [0.,0.01,0.03]+[wopT,]+[0.11,0.15]
# min , max are the same for both ...
wmin,wmax = min(wcf_plt_M),max(wcf_plt_M)
color = lambda www: plt.cm.brg(np.sqrt((www-wmin)/(wmax-wmin)))
colm = ["#843A36","#F0924F","#70831C","#046ABE","#A61E84","#E45A79"]
colt = ["#843A36","#F0924F","#70831C","#F71231","#A61E84","#E45A79"]
# plotting M, i.e. for bottom:
for i,wm in enumerate(wcf_plt_M):
    if wm not in data_sorted[(data_sorted.Model==mod)&(data_sorted.FFTemp==fft)]['W_coeff'].unique():
        pass
    else:
        #start plotting ...
        TM_data_wcf = cors_TM[cors_TM['W_coeff']==wm]
        x, y = TM_data_wcf['Temp'], TM_data_wcf['exp_M']
        ax_bottom.plot(x,y,color=colm[i],marker='o',ms=6,lw=2,mew=0.0,label="%.2f"%wm)
for i,wt in enumerate(wcf_plt_T):
    if wt not in data_sorted[(data_sorted.Model==mod)&(data_sorted.FFTemp==fft)]['W_coeff'].unique():
        pass
    else:
        #start plotting ...
        TM_data_wcf = cors_TM[cors_TM['W_coeff']==wt]
        x, y = TM_data_wcf['Temp'], TM_data_wcf['exp_T']
        ax_top.plot(x,y,color=colt[i],marker='o',ms=6,lw=2,mew=0.0,label="%.2f"%wt)
#
ax_top.set_xlim((0.35,2.1))
ax_bottom.set_xlim((0.35,2.1))
#
ax_top.yaxis.set_ticks_position('left')
ax_top.xaxis.set_ticks_position('bottom')
ax_bottom.yaxis.set_ticks_position('left')
ax_bottom.xaxis.set_ticks_position('bottom')
#
leg_top=ax_top.legend(loc='upper center',bbox_to_anchor=(0.9,1.0),numpoints=1,frameon=False,handlelength=0.8,handletextpad=0.25,labelspacing=0.75,borderaxespad=0.25,title=r"$w$")
leg_bottom=ax_bottom.legend(loc='upper center',bbox_to_anchor=(0.9,1.0),numpoints=1,frameon=False,handlelength=0.8,handletextpad=0.25,labelspacing=0.75,borderaxespad=0.25,title=r"$w$")
#
top_title = leg_top.get_title()
bottom_title = leg_bottom.get_title()
top_title.set_fontsize(11)
top_title.set_fontweight('bold')
bottom_title.set_fontsize(11)
bottom_title.set_fontweight('bold')
#
ax_top.xaxis.set_tick_params(labeltop='off')
ax_top.xaxis.set_tick_params(labelbottom='off')
#
ax_bottom.set_xlabel(r'$T$, p.u.',labelpad = 2)
ax_bottom.set_ylabel(r'$R_M$',rotation='horizontal',labelpad = 9)
ax_top.set_ylabel(r'$R_T$',rotation='horizontal',labelpad = 9)
#
ax_bottom.yaxis.label.set_size(9)
ax_top.yaxis.label.set_size(9)
#
ax_bottom.yaxis.set_major_locator( MaxNLocator(nbins = 7) )
ax_top.yaxis.set_major_locator( MaxNLocator(nbins = 7) )
#
# ax_top.vlines(1.0,-1,1)
ax_top.add_artist(ConnectionPatch(xyA=[T_T,RMT_max.RT_max.max()], xyB=[T_T,ax_top.get_ylim()[0]], coordsA="data", coordsB="data", axesA=ax_top, axesB=ax_top,color="#F71231",lw=2,zorder=100,ls='dashed'))
ax_bottom.add_artist(ConnectionPatch(xyA=[T_T,ax_bottom.get_ylim()[1]], xyB=[T_T,ax_bottom.get_ylim()[0]], coordsA="data", coordsB="data", axesA=ax_bottom, axesB=ax_bottom,color="#F71231",lw=2,zorder=100,ls='dashed'))
ax_bottom.add_artist(ConnectionPatch(xyA=[T_M,RMT_max.RT_max.max()], xyB=[T_M,ax_bottom.get_ylim()[0]], coordsA="data", coordsB="data", axesA=ax_bottom, axesB=ax_bottom,color="#046ABE",lw=2,zorder=100,ls='dashed'))
#
ax_bottom.text(T_M-0.1*(T_T-T_M),ax_bottom.get_ylim()[0]+0.00,r'$T_M$',fontsize=11,verticalalignment='bottom',horizontalalignment='right',color="#046ABE",fontweight='bold')
ax_bottom.text(T_T-0.1*(T_T-T_M),ax_bottom.get_ylim()[0]+0.00,r'$T_T$',fontsize=11,verticalalignment='bottom',horizontalalignment='right',color="#F71231",fontweight='bold')
#
width = 0.16
rectangle = Rectangle((0.9-0.5*width, 0.54), width, 0.075, fill=False, lw=2,color="#F71231", transform=ax_top.transAxes, zorder=20)
rectangle2 = Rectangle((0.9-0.5*width, 0.54), width, 0.075, fill=False, lw=2,color="#046ABE", transform=ax_bottom.transAxes, zorder=20)
ax_top.add_patch(rectangle)
ax_bottom.add_patch(rectangle2)
#
#
###################################################################################################
#  A & D FIGURE ...
#################################
# bottom axes 
ax_down = ax_TMAD[1,1]
# top axes 
ax_up = ax_TMAD[0,1]
################################
################################
################################
wcf_plt_A = [0.,0.01,0.03]+[wopMT,]+[0.11,0.15]
wcf_plt_D = [0.,0.01,0.03]+[wopMT,]+[0.11,0.15]
# min , max are the same for both ...
wmin,wmax = min(wcf_plt_A),max(wcf_plt_A)
color = lambda www: plt.cm.brg(np.sqrt((www-wmin)/(wmax-wmin)))
cola = ["#843A36","#F0924F","#70831C","#046ABE","#A61E84","#E45A79"]
# cold = ["#843A36","#F0924F","#70831C","#F71231","#A61E84","#E45A79"]
cold = ["#843A36","#F0924F","#70831C","#046ABE","#A61E84","#E45A79"]
# plotting M, i.e. for bottom:
for i,wa in enumerate(wcf_plt_A):
    if wa not in data_sorted[(data_sorted.Model==mod)&(data_sorted.FFTemp==fft)]['W_coeff'].unique():
        pass
    else:
        #start plotting ...
        TM_data_wcf = cors_TM[cors_TM['W_coeff']==wa]
        x, y = TM_data_wcf['Temp'], TM_data_wcf['exp_A']
        ax_up.plot(x,y,color=cola[i],marker='o',ms=6,lw=2,mew=0.0,label="%.2f"%wa)
for i,wt in enumerate(wcf_plt_D):
    if wt not in data_sorted[(data_sorted.Model==mod)&(data_sorted.FFTemp==fft)]['W_coeff'].unique():
        pass
    else:
        #start plotting ...
        TM_data_wcf = cors_TM[cors_TM['W_coeff']==wt]
        x, y = TM_data_wcf['Temp'], TM_data_wcf['exp_D']
        ax_down.plot(x,y,color=cold[i],marker='o',ms=6,lw=2,mew=0.0,label="%.2f"%wt)
#
ax_down.set_xlim((0.35,2.1))
ax_up.set_xlim((0.35,2.1))
#
ax_down.yaxis.set_ticks_position('left')
ax_down.xaxis.set_ticks_position('bottom')
ax_up.yaxis.set_ticks_position('left')
ax_up.xaxis.set_ticks_position('bottom')
#
#
leg_top=ax_down.legend(loc='upper center',bbox_to_anchor=(0.9,1.0),numpoints=1,frameon=False,handlelength=0.8,handletextpad=0.25,labelspacing=0.75,borderaxespad=0.25,title=r"$w$")
leg_bottom=ax_up.legend(loc='upper center',bbox_to_anchor=(0.9,1.0),numpoints=1,frameon=False,handlelength=0.8,handletextpad=0.25,labelspacing=0.75,borderaxespad=0.25,title=r"$w$")
#
top_title = leg_top.get_title()
bottom_title = leg_bottom.get_title()
top_title.set_fontsize(11)
top_title.set_fontweight('bold')
bottom_title.set_fontsize(11)
bottom_title.set_fontweight('bold')
#
ax_up.xaxis.set_tick_params(labeltop='off')
ax_up.xaxis.set_tick_params(labelbottom='off')
#
ax_up.set_ylabel('$R_A$',rotation='horizontal')#,labelpad = 9)
ax_down.set_ylabel('$R_D$',rotation='horizontal')#,labelpad = 5)
#
ax_up.yaxis.set_label_coords(-0.12, 0.5, transform=ax_up.transAxes)
ax_down.yaxis.set_label_coords(-0.12, 0.5, transform=ax_down.transAxes)
#
ax_down.set_xlabel('$T$, p.u.',labelpad = 2)
#
ax_up.yaxis.label.set_size(9)
ax_down.yaxis.label.set_size(9)
#
ax_up.yaxis.set_major_locator( MaxNLocator(nbins = 7) )
ax_down.yaxis.set_major_locator( MaxNLocator(nbins = 7) )
# #
# #
# # the x coords of this transformation are data, and the
# # y coord are axes
# trans_up = transforms.blended_transform_factory(ax_up.transData, ax_up.transAxes)
# trans_down = transforms.blended_transform_factory(ax_down.transData, ax_down.transAxes)
# # We want x to be in data coordinates and y to
# # span from 0..1 in axes coords
# rect_up = patches.Rectangle((Top_M,0), width=(Top_T-Top_M), height=1, transform=trans_up, color='yellow', alpha=0.5)
# rect_down = patches.Rectangle((Top_M,0), width=(Top_T-Top_M), height=1, transform=trans_down, color='yellow', alpha=0.5)
# rect_down_manual = patches.Rectangle((1.2,0), width=(1.5-1.2), height=1, transform=trans_down, color='blue', alpha=0.5)
# #
# #
ax_up.axvspan(Top_M, Top_T, facecolor='yellow', alpha=0.5)
ax_down.axvspan(Top_M, Top_T, facecolor='yellow', alpha=0.5)
#
plt.savefig(os.path.join(results_path,'%s_Figure_4.png'%exp_fname.replace('.','_')),dpi=300)










###########################################################
#  now, figure 4, or Rmax(w) for T,M and T+M
###########################################################
# RMT_max = pd.DataFrame(RMT_max,columns=['W_coeff','RM_max','RT_max'])
plt.clf()
x_fig_size = 3.34
v_coeff = 1.0
fig = plt.figure(figsize=(x_fig_size,v_coeff*x_fig_size))
# axes info ...
left = 0.12
bottom = 0.12
width = 0.95 - left
height = 0.97 - bottom
# bottom axes 
ax = plt.axes([left, bottom, width, height])

# # combination ...
# # RMT_max['RMT_max'] = RMT_max['RT_max']+RMT_max['RM_max']
# RMT_max['RMT_max'] = RMT_max['RT_max']/RMT_max['RT_max'].max()+RMT_max['RM_max']/RMT_max['RM_max'].max()
# wopMT = RMT_max['W_coeff'][RMT_max['RMT_max'].argmax()]
# ##############
ax.plot(RMT_max['W_coeff'],RMT_max['RM_max'],color='blue',marker='o',ms=6,lw=2,mew=0.0,label="$R_M$")
ax.plot(RMT_max['W_coeff'],RMT_max['RT_max'],color='red',marker='o',alpha=1.0,ms=6,lw=2,mew=0.0,label="$R_T$")
ax.plot(RMT_max['W_coeff'],RMT_max['RMT_max'],color='black',marker='o',ms=6,lw=2,mew=0.0,label="$R_M+R_T$")


ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')


ax.legend(loc='best',numpoints=1,frameon=False,handlelength=0.8,handletextpad=0.25,labelspacing=0.75,borderaxespad=0.25)

ax.xaxis.set_tick_params(labeltop='off')

ax.set_xlabel('$w$, adjustment parameter',labelpad = 2)
ax.set_ylabel('$R$, correlation--based metric',rotation='vertical',labelpad = 2)

ax_bottom.yaxis.label.set_size(9)
ax_top.yaxis.label.set_size(9)


ax.add_artist(ConnectionPatch(xyA=[wopT,RMT_max.RT_max.max()], xyB=[wopT,ax.get_ylim()[0]], coordsA="data", coordsB="data", axesA=ax, axesB=ax,color="red",lw=2,zorder=100,ls='dashed'))
ax.add_artist(ConnectionPatch(xyA=[wopM,RMT_max.RM_max.max()], xyB=[wopM,ax.get_ylim()[0]], coordsA="data", coordsB="data", axesA=ax, axesB=ax,color="blue",lw=2,zorder=100,ls='dashed'))
ax.add_artist(ConnectionPatch(xyA=[wopMT,RMT_max['RMT_max'].max()], xyB=[wopMT,RMT_max.RM_max.max()], coordsA="data", coordsB="data", axesA=ax, axesB=ax,color="black",lw=2,zorder=100,ls='dashed'))


ax.text(wopM-0.002,ax.get_ylim()[0]+0.02,'$w_{M}$',fontsize=11,verticalalignment='bottom',horizontalalignment='right',color="blue")
ax.text(wopT+0.002,ax.get_ylim()[0]+0.02,'$w_{T}$',fontsize=11,verticalalignment='bottom',horizontalalignment='left',color="red")
ax.text(wopMT-0.002,RMT_max.RM_max.max()+0.02,'$w^{*}$',fontsize=11,verticalalignment='bottom',horizontalalignment='right',color="black")


# ax_bottom.yaxis.set_major_locator( MaxNLocator(nbins = 7) )
# ax_top.yaxis.set_major_locator( MaxNLocator(nbins = 7) )

# ax_top.add_artist(ConnectionPatch(xyA=[T_T,RMT_max.RT_max.max()], xyB=[T_T,ax_bottom.get_ylim()[0]], coordsA="data", coordsB="data", axesA=ax_top, axesB=ax_bottom,color="#F71231",lw=2,zorder=100,ls='dashed'))
# ax_bottom.add_artist(ConnectionPatch(xyA=[T_M,RMT_max.RT_max.max()], xyB=[T_M,ax_bottom.get_ylim()[0]], coordsA="data", coordsB="data", axesA=ax_bottom, axesB=ax_bottom,color="#046ABE",lw=2,zorder=100,ls='dashed'))

# ax_bottom.text(T_M,ax_bottom.get_ylim()[0],'$T_M$',fontsize=11,verticalalignment='bottom',horizontalalignment='right',color="#046ABE")
# ax_bottom.text(T_T,ax_bottom.get_ylim()[0],'$T_T$',fontsize=11,verticalalignment='bottom',horizontalalignment='right',color="#F71231")

plt.savefig(os.path.join(results_path,'%s_SuppFigure4.png'%exp_fname),dpi=300)





# p = plt.axvspan(1.25, 1.55, facecolor='g', alpha=0.5)


#################################
# the brooms and their insets ...
plt.clf()
x_fig_size = 7.4
v_coeff = 0.45
fig = plt.figure(figsize=(x_fig_size,v_coeff*x_fig_size))
# between axes ...
hor_space = 0.02
# axes info ...
left = 0.06
bottom = 0.12
width = 0.5*(0.9 - left - hor_space)
height = 0.975 - bottom
# bottom axes 
ax_left = plt.axes([left, bottom, width, height])
left += (width + hor_space)
# top axes 
ax_right = plt.axes([left, bottom, width, height])
#############################################################
#############################################################
#############################################################
#############################################################
w_left = 0.0
#############################################################
left_dat = data_sorted[(data_sorted['W_coeff']==w_left)&(data_sorted.Model == 64)&(data_sorted.FFTemp == 'MJ99')]
for aa,col in zip(aacids,cols):
    ax_left.plot(left_dat['Temp'],left_dat[aa]*100.0,'o-',ms=4, lw=2, mew=0.0, color=col, label=aa)
#
# ax_left.yaxis.set_ticks_position('left')
ax_left.yaxis.set_ticks_position('both')
ax_left.xaxis.set_ticks_position('bottom')
#
ax_left.set_xlabel(r'T, p.u.')
ax_left.set_ylabel(r'$f_{a}$, \%')
#############################################################
w_right = 0.05
right_dat = data_sorted[(data_sorted['W_coeff']==w_right)&(data_sorted.Model == 64)&(data_sorted.FFTemp == 'MJ99')]
for aa,col in zip(aacids,cols):
    ax_right.plot(right_dat['Temp'],right_dat[aa]*100.0,'o-',ms=4, lw=2, mew=0.0, color=col, label=r'\texttt{%s}'%aa)
################################
ax_right.yaxis.set_ticks_position('left')
ax_right.xaxis.set_ticks_position('bottom')
################################
ax_right.set_xlabel(r'T, p.u.')
# ax_right.set_ylabel('$f_{a}$, \%')
################################
leg = ax_right.legend(loc='upper left',numpoints=1, frameon=False, bbox_to_anchor=(1.01, 1.035),fontsize=8,handlelength=1.0,markerscale=1.5,labelspacing=0.35)
# ax.legend(loc='best',numpoints=1,frameon=False,handlelength=0.8,handletextpad=0.25,labelspacing=0.75,borderaxespad=0.25)
#
# ax_right.axvspan()
#
ax_right.axvspan(Top_M, Top_T, facecolor='yellow', alpha=0.5)
#
# # leg = plt.legend(loc='best')
# # leg.get_frame().set_alpha(0)
# # leg.get_frame().set_edgecolor('white')
# # wider markers in the legend ...
# for legobj in leg.legendHandles:
#     legobj.set_linewidth(3.0)
###########################################
###########################################
###########################################
ax_right.yaxis.set_tick_params(labelleft='off')
ax_left.yaxis.set_tick_params(labelright='off')
# ax_right.xaxis.set_tick_params(labelbottom='off')
#
ax_left.set_ylim((0,18))
ax_right.set_ylim((0,18))
#
ax_left.set_xlim((0.35,1.75))
ax_right.set_xlim((0.35,1.75))
##########################################
# inset plotting ...
###########################################
left -= (width + hor_space)
inset_left = fig.add_axes([left+0.048, 0.65, 0.2, 0.3],zorder=500)
left += (width + hor_space)
inset_right = fig.add_axes([left+0.048, 0.65, 0.2, 0.3],zorder=500)


left_charged = left_dat[list('DEKR')].sum(axis=1)
left_hydrophob = left_dat[list('AGNQSTHY')].sum(axis=1)
left_hydrophil = left_dat[list('MPCLVWIF')].sum(axis=1)    

right_charged = right_dat[list('DEKR')].sum(axis=1)
right_hydrophob = right_dat[list('AGNQSTHY')].sum(axis=1)
right_hydrophil = right_dat[list('MPCLVWIF')].sum(axis=1)    


inset_left.plot(left_dat.Temp[::2],left_hydrophob[::2]*100.0,'o-',ms=6, lw=2, mew=0.0, color='g', label=r'\texttt{AGNQSTHY}')
inset_left.plot(left_dat.Temp[::2],left_hydrophil[::2]*100.0,'o-',ms=6, lw=2, mew=0.0, color='b', label=r'\texttt{MPCLVWIF}')
inset_left.plot(left_dat.Temp[::2],left_charged[::2]*100.0,'o-',ms=6, lw=2, mew=0.0, color='r', label=r'\texttt{DEKR}')


inset_right.plot(right_dat.Temp[::2],right_hydrophob[::2]*100.0,'o-',ms=6, lw=2, mew=0.0, color='g', label=r'\texttt{AGNQSTHY}')
inset_right.plot(right_dat.Temp[::2],right_hydrophil[::2]*100.0,'o-',ms=6, lw=2, mew=0.0, color='b', label=r'\texttt{MPCLVWIF}')
inset_right.plot(right_dat.Temp[::2],right_charged[::2]*100.0,'o-',ms=6, lw=2, mew=0.0, color='r', label=r'\texttt{DEKR}')


inset_left.yaxis.set_major_locator( MaxNLocator(nbins = 4) )
inset_left.xaxis.set_major_locator( MaxNLocator(nbins = 4) )
#
inset_right.yaxis.set_major_locator( MaxNLocator(nbins = 4) )
inset_right.xaxis.set_major_locator( MaxNLocator(nbins = 4) )


# inset_left.xaxis.set_tick_params(labeltop='off')
# inset_left.xaxis.set_tick_params(labelbottom='off')


inset_left.tick_params(axis='x',which='both',top='off',bottom='on',pad=2.5)
inset_left.tick_params(axis='y',which='both',left='on',right='off',pad=2.5)
#
inset_right.tick_params(axis='x',which='both',top='off',bottom='on',pad=2.5)
inset_right.tick_params(axis='y',which='both',left='on',right='off',pad=2.5)

inset_left.set_xlim((0.35,1.65))
inset_right.set_xlim((0.35,1.65))
#
inset_left.set_ylim((18,46))
inset_right.set_ylim((18,46))


inset_left.set_xlabel('T, p.u.',labelpad=1)
inset_right.set_xlabel('T, p.u.',labelpad=1)
#
lab_left = inset_left.set_ylabel(r'$f_{a}$, \%',labelpad=2)
lab_right = inset_right.set_ylabel(r'$f_{a}$, \%',labelpad=2)


lab_left.set_bbox(dict(facecolor='white', alpha=1.0, edgecolor='None'))
lab_right.set_bbox(dict(facecolor='white', alpha=1.0, edgecolor='None'))
# xxx.set_clip_on(True)

# leg = plt.legend(loc='center left',fontsize=15)
leg = inset_left.legend(loc='upper left',numpoints=1, bbox_to_anchor=(0.99, 1.02), frameon=False ,fontsize=10,handlelength=0.8,handletextpad=0.18,labelspacing=0.5,borderaxespad=0.1)
leg = inset_right.legend(loc='upper left',numpoints=1, bbox_to_anchor=(0.99, 1.02), frameon=False ,fontsize=10,handlelength=0.8,handletextpad=0.18,labelspacing=0.5,borderaxespad=0.1)
# # leg = plt.legend(loc='best')
# leg.get_frame().set_alpha(0)
# leg.get_frame().set_edgecolor('white')
# # wider markers in the legend ...
# for legobj in leg.legendHandles:
#     legobj.set_linewidth(6.0)

# # ax1.text(-0.15, 1.15, 'A', transform=ax1.transAxes, fontsize=10, weight='bold', va='top', color='black')
# # ax2.text(-0.15, 1.15, 'B', transform=ax2.transAxes, fontsize=10, weight='bold', va='top', color='black')
###########################################
###########################################
fig.savefig(os.path.join(results_path,"Figure_3.png"),dpi=300)
###########################################
###########################################







###############################################
###############################################
##### SLopes vs evolved slopes ...
###############################################
###############################################
M=64
FFT='MJ99'
W = wopMT
# M,W=MW_case
plt.clf()
###########
x_fig_size = 3.34
v_coeff = 0.85
fig, ax = plt.subplots(1,1,figsize=(x_fig_size,v_coeff*x_fig_size))#, sharex=True, sharey=True)
###########
sub_data = data_sorted[(data_sorted.Model == M)&(data_sorted.FFTemp == FFT)&(data_sorted.W_coeff == W)].copy()
# data used to be offsetted by the quantity of rows ....
# offset_left = 25 if sub_data.shape[0]>=29 else int(sub_data.shape[0]*0.55)
# offset_right = 2 if sub_data.shape[0]>=29 else int(sub_data.shape[0]*0.1)+1
# sub_off_data = sub_data[offset_left:-offset_right]
# not any more - visual fitting check is applied ...
sub_off_data = sub_data[(sub_data.Temp>=Top_M)&(sub_data.Temp<=Top_T)]
data_TT = sub_data[(sub_data.Temp==Top_T)]
data_TM = sub_data[(sub_data.Temp==Top_M)]
# sub_off_data = sub_data[(sub_data.Temp>=T_M)&(sub_data.Temp<=T_T)]
#
sim_D = sub_off_data[aacids].apply(lambda x: x.cov(sub_off_data.Temp)/sub_off_data.Temp.var())
sim_D2 = (data_TT[aacids].reset_index(drop=True) - data_TM[aacids].reset_index(drop=True))/(Top_T-Top_M)
sim_D2.convert_objects(convert_numeric=True)
sim_D = sim_D2.T[0]
# sim_D = sub_off_data[aacids].apply(lambda x: x.cov(sub_off_data.Temp)/sub_off_data.Temp.var())
#
for lab,xx,yy in zip(exp_D.index,exp_D,sim_D):
    ax.scatter(xx, yy,s=3, marker="o",c='white',edgecolors='none')
    ax.text(xx,yy,r"\texttt{%s}"%lab,color='black',horizontalalignment='center',verticalalignment='center',fontsize=13,fontweight='bold')
#
# ax.plot(exp_D,sim_D,'ro',label='slopes M%d FFT%s W%.4f'%(M,FFT,W))
a,b,r,pval,_ = st.linregress(exp_D.values,sim_D.values)
exp_D_range = np.linspace(exp_D.min()*1.01,exp_D.max()*1.01)
ax.plot(exp_D_range,a*exp_D_range+b,color='silver',linewidth=1.5,linestyle='-',zorder=100,label=r'linear fit: R=%.2f, $p=%.3f$'%(r,pval))
# ax.plot(exp_D_range,exp_D_range,color='red',linewidth=1.0,linestyle='--',zorder=102)
# ax.set_title("sim_D = a*exp_D+b, a=%.4f, b=%.4f, r=%.2f"%(a,b,r))
ax.legend(loc='upper left',frameon=False,handlelength=1.5,handletextpad=0.1)
# ax.legend(loc='upper left',bbox_to_anchor=(0.,0.9),frameon=False,handlelength=1.5,handletextpad=0.1)
# ax.legend(loc='best',frameon=False)
ax.set_xlabel(slopes_label)
ax.set_ylabel("simulated slopes, 1/p.u.")
# ax.set_xlabel(slopes_label,labelpad=1.5)
# ax.set_ylabel("simulated slopes, 1/p.u.",labelpad=0.5)
plt.tight_layout(pad=0.4, h_pad=None, w_pad=None)
#
# '$f_{a}$, \%'
# signif_cai10_bacter
# top 10% CAI t.o. bacteria
# top 10% CAI t.o. archaea
# proteome bacteria
# proteome archaea
#
ax.xaxis.set_major_locator( MaxNLocator(nbins = 6) )
ax.yaxis.set_major_locator( MaxNLocator(nbins = 6) )
#
ax.tick_params(axis='x',which='both',top='off',bottom='on',pad=3)
ax.tick_params(axis='y',which='both',left='on',right='off',pad=3)
#
# # ax.set_xlim((-0.025,0.03))
# # ax.set_ylim((-0.04,0.081))
ymin,ymax = sim_D.min(),sim_D.max()
y_span = ymax-ymin
ax.set_ylim((ymin-0.1*y_span,ymax+0.3*y_span))
# ymax,ymin = ax.get_ylim()
# xmax,xmin = ax.get_xlim()
# ax.plot([0,0],[-0.1*(ymax-ymin),0],color='blue',linewidth=1.0,linestyle='-',zorder=104)
# ax.plot([-0.1*(xmax-xmin),0],[0,0],color='blue',linewidth=1.0,linestyle='-',zorder=104)
# #
# plt.legend(loc='best')
fig.savefig(os.path.join(results_path,"%s_Slopes_sim_vs_exp_protbact.png"%exp_fname),dpi=300)













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
# up to now, we agreed to use Model=64 and FFT=MJ99 only, but anyways ...
mod = 64
fft = data_sorted[data_sorted.Model==mod].FFTemp.unique()[0]
####################################################
data_mod = data_sorted[(data_sorted.Model==mod)&(data_sorted.FFTemp==fft)]
wcf_desired = [0.,0.02,0.04,0.05,0.06,0.07,0.09,0.12,0.15]+[wopMT]
# get reid of the repeats if neccessary ...
wcf_desired = sorted(list(set(wcf_desired)))
grouped = data_mod.groupby('W_coeff')
for i,(wcf,dat) in enumerate(grouped):
    if wcf in wcf_desired:
        datt = dat.transpose()
        ddt = datt[:20] #amino acids only ...
        proteome_argentina_cost = argentina_cost.T.dot(ddt)
        if wcf == wopMT:
            ax_right.plot(datt.loc['Temp'],proteome_argentina_cost.T,color=cols[i],marker='s',ms=8,mew=0,lw=2,label="%.2f"%wcf)
        else:
            ax_right.plot(datt.loc['Temp'],proteome_argentina_cost.T,color=cols[i],marker='o',ms=6,mew=0,lw=2,label="%.2f"%wcf)
##################################
ax_left.yaxis.set_ticks_position('left')
ax_left.xaxis.set_ticks_position('bottom')
# leg = ax_left.legend(loc='best',frameon=False,numpoints=1,markerscale=2,fontsize=13)
# leg.get_frame().set_alpha(0)
# leg.get_frame().set_edgecolor('white')
# ax_left.set_xlim((0.0,6.1))
# plt.tight_layout()
# ax_left.set_title("Arg_proteome_cost_mod%d_fft%s.ghpcc.png"%(mod,fft))
# fig.savefig(os.path.join(results_path,"Arg_proteome_cost_mod%d_fft%s.ghpcc.png"%(mod,fft)))
#####################################
#####################################
#
data_mod = data_sorted[(data_sorted.Model==mod)&(data_sorted.FFTemp==fft)]
grouped = data_mod.groupby('W_coeff')
for i,(wcf,dat) in enumerate(grouped):
    if wcf in wcf_desired:
        datt = dat.transpose()
        ddt = datt[:20] #amino acids only ...
        proteome_akashi_cost = akashi_cost.T.dot(ddt)
        if wcf == wopMT:
            ax_left.plot(datt.loc['Temp'],proteome_akashi_cost.T,color=cols[i],marker='s',ms=8,mew=0,lw=2,label="%.2f"%wcf)
        else:
            ax_left.plot(datt.loc['Temp'],proteome_akashi_cost.T,color=cols[i],marker='o',ms=6,mew=0,lw=2,label="%.2f"%wcf)
##################################
ax_right.yaxis.set_ticks_position('left')
ax_right.xaxis.set_ticks_position('bottom')
# leg = ax_right.legend(loc='best',frameon=False,numpoints=1,markerscale=2,fontsize=13)
leg = ax_right.legend(loc='upper left',numpoints=1, frameon=False, fontsize=9, bbox_to_anchor=(1.01, 1.035),handlelength=1.5,markerscale=1.5,labelspacing=0.9, title=r"$w$")
#
#
ax_left.set_xlabel(r'T, p.u.',labelpad=1)
ax_right.set_xlabel(r'T, p.u.',labelpad=1)
#
#
ax_right.set_ylabel('AA maintenance cost, ATP/time',labelpad=1)
ax_left.set_ylabel('AA synthesis cost, ATP',labelpad=1)
#
#
ax_left.axvspan(Top_M, Top_T, facecolor='yellow', alpha=0.5)
ax_right.axvspan(Top_M, Top_T, facecolor='yellow', alpha=0.5)
# leg.get_frame().set_alpha(0)
# leg.get_frame().set_edgecolor('white')
# ax_right.set_xlim((0.0,6.1))
ax_left.set_xlim((0.4,1.7))
ax_right.set_xlim((0.4,1.7))
ax_right.set_ylim((25,230))
# plt.tight_layout()
# ax_right.set_title("Akashi_proteome_cost_mod%d_fft%s.ghpcc.png"%(mod,fft))
fig.savefig(os.path.join(results_path,"%s_Fig7_costs.png"%exp_fname),dpi=300)






























































