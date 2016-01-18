import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess as sub
import pandas as pd
import numpy as np
import sys
import re
import os
import scipy.stats as st

import random as rnd

from matplotlib.ticker import MaxNLocator
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Rectangle

import scipy.interpolate as interpol
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


# constants ...
# ALPHABET=20
# amino acis alphabet ...
# aacids = data.columns.values[:ALPHABET]
aacids = sorted('CMFILVWYAGTSNQDEHRKP')


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


# getting file name based on criteria combination ...
def get_exp_fname(cds_criteria,org_criteria,kingdom):
    if kingdom in ['arch','archaea']:
        return "exp_MTAD_%s_%s.arch.summary"%(cds_criteria,org_criteria)
    elif kingdom in ['bact','bacteria']:
        return "exp_MTAD_%s_%s.bact.summary"%(cds_criteria,org_criteria)
    else:
        raise TypeError('only archaeal and bacterial kingdoms are supported!')



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
data_exp = pd.read_csv(exp_fname,index_col=0)
# loaded ...
# give em names, as I figured, indexes and column names are awesome things ...
exp_T = data_exp['exp_T']
exp_M = data_exp['exp_M']
exp_A = data_exp['exp_A']
exp_D = data_exp['exp_D']
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
shuffle_fname = os.path.join(simul_path,"shuffled.dat")
shuffle_slopes_fname = os.path.join(simul_path,"shuffled_slopes.dat")
# loading simul data ...
data_sorted_D = pd.read_csv(shuffle_slopes_fname)
data_sorted_A = pd.read_csv(shuffle_fname)
# loaded ...
########################################################
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


######################################################################################
######################################################################################
######################################################################################
######################################################################################
plt.clf()
x_fig_size = 7.3
v_coeff = 0.45
fig = plt.figure(figsize=(x_fig_size,v_coeff*x_fig_size))
# between axes ...
hor_space = 0.07
# axes info ...
left = 0.06
bottom = 0.12
width = 0.5*(0.9 - left - hor_space)
height = 0.98 - bottom
# bottom axes 
ax_left = plt.axes([left, bottom, width, height])
left += (width + hor_space)
# top axes 
ax_right = plt.axes([left, bottom, width, height])
#
#
#
###############################################################
###############################################################
###############################################################
###############################################################
# shuff histo ...
t_regime = 'A'
mod = 64
# wcf = 0.02


data_sorted_D = data_sorted_D[data_sorted_D.Temp<1.6]
# #
# ##
# ###
# #####
# #######
# #######
rmin,rmax,num = -1.0,1.0,25
bins = np.linspace(rmin,rmax,num=num)
#########
for wcf,df_w in data_sorted_D.groupby('W_coeff'):
    corrs_w = []
    for shuf,df in df_w.groupby('Shuff'):
        sim_D = df[aacids].apply(lambda x: x.cov(df.Temp)/df.Temp.var())
        corrs_w.append(exp_D.corr(sim_D))
    if len(corrs_w)!= 100:
        print "something went wrong: number of shuffles isn't 100!"
        sys.exit(1)
    # get p value ...
    false_cases_num = sum(1 for cw in corrs_w[1:] if cw>corrs_w[0])
    pvalue = false_cases_num*1.0/len(corrs_w)
    #
    ################
    # do plotting ...
    ################
    if wcf == 0.06:
        print "plotting right plot: slopes"
        ax_right.hist(corrs_w, bins=bins,ec='none',normed=False, label='shuffled')
        width = (rmax-rmin)/(num-1)
        wt_bin_num = int((corrs_w[0]-rmin)/width)
        ax_right.bar([wt_bin_num*width+rmin,],[1.0,],width=width,color='red',edgecolor='None',label='predicted')
        # ax_right.bar([corrs_w[0]-0.5*width,],[1.0,],width=width,color='red',edgecolor='None',label='predicted')
        ax_right.set_xlim((rmin,rmax))
        ax_right.yaxis.set_ticks_position('left')
        ax_right.xaxis.set_ticks_position('bottom')
        ax_right.set_xlabel(r"$R_D$, slopes correlation")
        f_title = lambda pval: (r"p=%.3f"%pval) if (pval>=0.001) else r"p\textless0.001"
        leg_right = ax_right.legend(loc='upper right',frameon=False, title=f_title(pvalue))
        # ax_right.set_title('$w$=%.2f'%wcf)
        # ax_right.text(0.1,0.87,r'$w=%.2f$'%wcf,transform=ax_right.transAxes)
        #
        for legend_item in leg_right.get_patches():
            legend_item.set_edgecolor('none')
        # plt.title("Hist_mod%d_wcf%.3f_T%.1f.new_shuff2.pdf"%(mod,wcf,3.00))
        # fig.savefig(os.path.join(results_path,"Hist_%s_mod%d_wcf%.3f.new_shuff2.pdf"%(t_regime,mod,wcf)))
    else:
        pass




for wcf,df in data_sorted_A.groupby('W_coeff'):
    # shuff_at_wcf = data_sorted[data_sorted.W_coeff==wcf]
    data_shuff_T = df[aacids].transpose()
    # data_shuff_T_AA = data_shuff_T[:20]
    cors = data_shuff_T.convert_objects(convert_numeric=True).corrwith(exp_A)
    ##################
    # get p value ...
    false_cases_num = sum(1 for cw in cors.values[1:] if cw>cors.values[0])
    pvalue = false_cases_num*1.0/len(cors.values)
    #
    ##################
    # print wcf
    if (cors.shape[0]>1)and(wcf==0.06):
        print "plotting left plot: composition"
        ax_left.hist(cors.values, bins=bins,ec='none',normed=False,label='shuffled')
        ax_left.set_xlim((rmin,rmax))
        ax_left.yaxis.set_ticks_position('left')
        ax_left.xaxis.set_ticks_position('bottom')
        #
        width = (rmax-rmin)/(num-1)
        wt_bin_num = int((cors.values[0]-rmin)/width)
        ax_left.bar([wt_bin_num*width+rmin,],[1.0,],width=width,color='red',edgecolor='None',label='predicted')
        f_title = lambda pval: (r"p=%.3f"%pval) if (pval>=0.001) else r"p\textless0.001"
        leg_left = ax_left.legend(loc='upper right',frameon=False,title=f_title(pvalue))
        ax_left.set_xlabel(r"$R_A$, composition correlation")
        # ax_left.text(0.1,0.87,r'$w=%.2f$'%wcf,transform=ax_left.transAxes)
        #
        # set the same y limits ...
        counts_max_right = ax_right.get_ylim()[1]
        counts_max_left = ax_left.get_ylim()[1]
        ax_left.set_ylim( ( 0,max(counts_max_right,counts_max_left) ) )
        ax_right.set_ylim( ( 0,max(counts_max_right,counts_max_left) ) )
        # #
        #
        for legend_item in leg_left.get_patches():
            legend_item.set_edgecolor('none')
        #
        #
        ax_left.set_ylabel('histogram counts')
        #
        #
        ax_right.yaxis.set_tick_params(labelleft='off')
        ax_left.yaxis.set_tick_params(labelright='off')
        #
        #
        #
        # ax_left.text(0.05,0.87,'A',fontweight='bold',fontsize=18,transform=ax_left.transAxes)
        # ax_right.text(0.05,0.87,'B',fontweight='bold',fontsize=18,transform=ax_right.transAxes)
        #
        #
        #
        ax_left.tick_params(axis='x',which='both',top='off',bottom='on',pad=3)
        ax_left.tick_params(axis='y',which='both',left='on',right='off',pad=3)
        #
        ax_right.tick_params(axis='x',which='both',top='off',bottom='on',pad=3)
        ax_right.tick_params(axis='y',which='both',left='on',right='off',pad=3)
        #
        # # plt.legend(loc='best')
        # # plt.title("Hist_mod%d_wcf%.3f_T%.1f.new_shuff2.pdf"%(mod,wcf,3.00))
        # fig.savefig(os.path.join(results_path,"T_Hist_%s_mod%d_wcf%.3f.new_shuff2.pdf"%(t_regime,mod,wcf)))
    else:
        pass


fig.savefig(os.path.join(results_path,'%s_MainFig6.png'%exp_fname),dpi=300)


##############################
ra_wcf_df = {}
for wcf,df in data_sorted_A.groupby('W_coeff'):
    data_shuff_T = df[aacids].transpose()
    cors = data_shuff_T.convert_objects(convert_numeric=True).corrwith(exp_A)
    ##################
    cors.reset_index(drop=True,inplace=True)
    ra_wcf_df[wcf] = cors
######################
ra_wcf_df = pd.DataFrame(ra_wcf_df)
ra_wcf_df = ra_wcf_df.transpose()
#######################
#######################
plt.clf()
x_fig_size = 7.3
v_coeff = 0.45
fig = plt.figure(figsize=(x_fig_size,v_coeff*x_fig_size))
# between axes ...
hor_space = 0.09
# axes info ...
left = 0.1
bottom = 0.12
width = 0.5*(0.98 - left - hor_space)
height = 0.98 - bottom
# bottom axes 
ax_left = plt.axes([left, bottom, width, height])
left += (width + hor_space)
# top axes 
ax_right = plt.axes([left, bottom, width, height])
#######################
#######################
########
####
col_ff_log = lambda x,xmin,xmax: pd.np.log(x/xmin)/pd.np.log(xmax/xmin) 
col_ff_lin = lambda x,xmin,xmax: (x-xmin)/(xmax-xmin)+0.1
# col_ff_sq = lambda x,xmin,xmax: (x**2-xmin**2)/(xmax**2-xmin**2) 
col_ff_sq = lambda x,xmin,xmax: pd.np.sqrt(x-xmin)/pd.np.sqrt(xmax-xmin)
col_ff_pow = lambda x,xmin,xmax: pd.np.power((x-xmin),0.35)/pd.np.power((xmax-xmin),0.3)
col_norm = mpl.colors.Normalize(vmin=-0.03,vmax=0.12)
##
# left panel,
ramin,ramax = (lambda x: (x.max().min(),x.max().max()))(ra_wcf_df)
shuffs = range(100)
rnd.shuffle(shuffs)
for shuff in shuffs:
    ra_max = ra_wcf_df[shuff].max()
    wcf_max = ra_wcf_df[shuff].argmax()
    ax_left.plot(ra_wcf_df.index,ra_wcf_df[shuff],'o-',color=plt.cm.hot_r(col_norm(wcf_max)),mew=0)
#
#
ax_left.set_ylim((-1,1))
ax_left.set_xlim((-0.005,0.125))
ax_left.set_ylabel(r'$R_A$, correlation coefficient')
ax_left.set_xlabel(r'$w$, cost adjustment parameter')
#
#
#
# right panel,
# w_range = (lambda x: (x.min(),x.max()) )(ra_wcf_df.index)
# bins = pd.np.linspace(w_range[0],w_range[1],ra_wcf_df.index.size)
bins = pd.np.arange(0,0.14,0.02)-0.01
#
ax_right.hist(ra_wcf_df.apply(lambda x: x.argmax()),bins=bins,edgecolor='deepskyblue',color='dodgerblue')
ax_right.set_xlim((bins[0],bins[-1]))
ax_right.set_ylabel('histogram counts')
ax_right.set_xlabel(r'$w^*$, optimal parameter')

ax_right.set_ylim((0,43))
ax_right.yaxis.set_major_locator( MaxNLocator(nbins = 5) )

cax = fig.add_axes([0.11,0.26,0.3,0.04])
cmap = mpl.cm.hot_r
# norm = mpl.colors.Normalize(vmin=5, vmax=10)
cbar = mpl.colorbar.ColorbarBase(cax, cmap=mpl.cm.hot_r, norm=col_norm, boundaries = [0,0.02,0.04,0.06,0.08,0.1,0.12,0.14],orientation='horizontal')
cbar.set_ticks(pd.np.asarray([0.0,0.02,0.04,0.06,0.08,0.1,0.12])+0.01)
cbar.set_ticklabels(['0.0','0.02','0.04','0.06','0.08','0.1','0.12'])
cbar.set_label(r'$w^*$, optimal parameter',labelpad=1)
# cbar.set_clim(vmin=0.02,vmax=0.12)
# ax_right.yaxis.set_tick_params(labelleft='off')
# ax_left.yaxis.set_tick_params(labelright='off')
# cax.xaxis.set_tick_params
cax.tick_params(axis='x',which='both',top='on',bottom='on',length=2.5,labelsize=7.5)

ax_left.tick_params(axis='x',which='both',top='off',bottom='on')
ax_left.tick_params(axis='y',which='both',right='off',left='on')
#
ax_right.tick_params(axis='x',which='both',top='off',bottom='on')
ax_right.tick_params(axis='y',which='both',right='off',left='on')


# ax_left.text(0.031,0.925,'A',transform=ax_left.transAxes,fontsize=14,fontweight='bold')
# ax_right.text(0.031,0.925,'B',transform=ax_right.transAxes,fontsize=14,fontweight='bold')


fig.savefig(os.path.join(results_path,'%s_SuppFig3.png'%exp_fname),dpi=300)



# ######################################################################################
# ######################################################################################
# ######################################################################################
# ######################################################################################
# # SHUFFLED RESULTS FOR SUPPLEMENT ...
# ######################################################################################
# ######################################################################################
# ######################################################################################
# ######################################################################################
# plt.clf()
# x_fig_size = 7.3
# v_coeff = 0.45
# fig = plt.figure(figsize=(x_fig_size,v_coeff*x_fig_size))
# # between axes ...
# hor_space = 0.07
# # axes info ...
# left = 0.06
# bottom = 0.12
# width = 0.5*(0.9 - left - hor_space)
# height = 0.98 - bottom
# # bottom axes 
# ax_left = plt.axes([left, bottom, width, height])
# left += (width + hor_space)
# # top axes 
# ax_right = plt.axes([left, bottom, width, height])
# #
# #
# #
# ###############################################################
# ###############################################################
# ###############################################################
# ###############################################################
# # shuff histo ...
# t_regime = 'A'
# mod = 64
# #
# #
# #
# #
# # wcf = 0.02
# #
# #
# #
# #
# for wcf,df in data_sorted_A.groupby('W_coeff'):
#     # shuff_at_wcf = data_sorted[data_sorted.W_coeff==wcf]
#     data_shuff_T = df[aacids].transpose()
#     # data_shuff_T_AA = data_shuff_T[:20]
#     cors = data_shuff_T.convert_objects(convert_numeric=True).corrwith(exp_T)
#     ##################
#     ##################
#     # print wcf
#     #
#     #

#     #
#     #
#     #
#     #
#     if (cors.shape[0]>1)and(wcf==0.06):
#         print "plotting left plot: composition"
#         ax_left.hist(cors.values, bins=bins,ec='none',normed=False,label='shuffled')
#         ax_left.set_xlim((rmin,rmax))
#         ax_left.yaxis.set_ticks_position('left')
#         ax_left.xaxis.set_ticks_position('bottom')
#         #
#         width = (rmax-rmin)/(num-1)
#         wt_bin_num = int((cors.values[0]-rmin)/width)
#         ax_left.bar([wt_bin_num*width+rmin,],[1.0,],width=width,color='red',edgecolor='None',label='predicted')
#         leg_left = ax_left.legend(loc='best',frameon=False)
#         ax_left.set_xlabel("$R_{A}$, composition correlation")
#         ax_left.text(0.13,0.87,'($w$=%.2f)'%wcf,transform=ax_left.transAxes)
#         #
#         # set the same y limits ...
#         counts_max_right = ax_right.get_ylim()[1]
#         counts_max_left = ax_left.get_ylim()[1]
#         ax_left.set_ylim( ( 0,max(counts_max_right,counts_max_left) ) )
#         ax_right.set_ylim( ( 0,max(counts_max_right,counts_max_left) ) )
#         # #
#         #
#         for legend_item in leg_left.get_patches():
#             legend_item.set_edgecolor('none')
#         #
#         #
#         ax_left.set_ylabel('histogram counts')
#         #
#         #
#         ax_right.yaxis.set_tick_params(labelleft='off')
#         ax_left.yaxis.set_tick_params(labelright='off')
#         #
#         #
#         #
#         ax_left.text(0.05,0.87,'A',fontweight='bold',fontsize=18,transform=ax_left.transAxes)
#         ax_right.text(0.05,0.87,'B',fontweight='bold',fontsize=18,transform=ax_right.transAxes)
#         #
#         #
#         #
#         ax_left.tick_params(axis='x',which='both',top='off',bottom='on',pad=3)
#         ax_left.tick_params(axis='y',which='both',left='on',right='off',pad=3)
#         #
#         ax_right.tick_params(axis='x',which='both',top='off',bottom='on',pad=3)
#         ax_right.tick_params(axis='y',which='both',left='on',right='off',pad=3)
#         #
#         # # plt.legend(loc='best')
#         # # plt.title("Hist_mod%d_wcf%.3f_T%.1f.new_shuff2.pdf"%(mod,wcf,3.00))
#         # fig.savefig(os.path.join(results_path,"T_Hist_%s_mod%d_wcf%.3f.new_shuff2.pdf"%(t_regime,mod,wcf)))
#     else:
#         pass


# fig.savefig(os.path.join(results_path,'validation.pdf'))







