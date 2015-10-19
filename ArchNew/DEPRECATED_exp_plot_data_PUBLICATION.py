import os
import sys
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats as st
import random as rnd
#
#
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullFormatter
import scipy.interpolate as interpol
font = {'family' : 'sans-serif',
        #'weight' : 'bold',
        'size'   :9}
#
mpl.rc('font', **font)
#
#
#
#
#
#
#
#####
########
############
###############
# akashi-cost.d
# argentina-cost.d
# "/Users/venevs/Desktop/Dropbox/protein_design/NEW_DATASETS"
# mac_root_path = r"/Users/venevs/Dropbox (Personal)"
# macbook_root_path = r"/Users/sergpolly/Desktop/Dropbox (Personal)"
# cosmos_root_path = "/home/venevs/Dropbox"
# /Users/venevs/Dropbox (Personal)/protein_design/Correct_MJ99_Argentina_PUBLICATION
# root_path = macbook_root_path
root_path = "."
# genome_data_path = os.path.join(root_path,"protein_design/NEW_DATASETS/ftp.ncbi.nih.gov/genomes/all")
# data_file = os.path.join(root_path,"protein_design","NEW_DATASETS/Sp_Temp_Compos_Phylog.dat")
# datra ready to use ...
# data_file = "./Sp_Temp_Compos_Phylog.dat"
data_file = "./arch_whole_proteome.dat"
# data_file_bact = "./EUBACTER_83.dat"

results_path = "."
# os.mkdir(results_path)
dat = pd.read_csv( data_file )
# dat_bact = pd.read_csv( data_file_bact )
##############################
aacids = ['C', 'M', 'F', 'I', 'L', 'V', 'W', 'Y', 'A', 'G', 'T', 'S', 'N', 'Q', 'D', 'E', 'H', 'R', 'K', 'P']
# cm = plt.get_cmap('gist_rainbow')
# different slices of data ...
# indcs = (dat.Q!=0.0)&(dat.GC>=30)&(dat.GC<=65)
# # indcs = (dat.aa_Q!=0.0)&(dat.subdivision!='Halobacteria')
# indcs = (dat.Q!=0.0)&(dat.subdivision!='Halobacteria')
# indcs_all_reasonable = indcs
# indcs_methan = (dat.Q!=0.0)&(dat.subdivision!='Halobacteria')&(dat.subdivision=='Methanococci')
# indcs_GC50 = (dat.Q!=0.0)&(dat.subdivision!='Halobacteria')&(dat.GC>=45)&(dat.GC<=55)
# indcs_GC30 = (dat.Q!=0.0)&(dat.subdivision!='Halobacteria')&(dat.GC>=26)&(dat.GC<=34)
######
# #
# # bacteria by Goncearenco ...
# # let's undersample parasites a bit, manually ...
# not_parasites = list((dat_bact.topt != 37).nonzero()[0])
# parasites = list((dat_bact.topt == 37).nonzero()[0])
# parasites8 = list(rnd.sample(parasites,6))
# #
# dat_bact = dat_bact.iloc[not_parasites + parasites8]
# dat_bact[aacids] = dat_bact[aacids]*100.0
# #
# #
# # use only interesting data ...
dat[aacids] = dat[aacids]*100.0
######
GCmin,GCmax = 20,70
tmin,tmax = 15,110
aamin,aamax = 0,12
######
#########################
#########################
#########################
# # prokaryotes 250 has no GC information, but might interesting for the final table (imshow thing) ...
# prok250_old_KZ_fname = os.path.join(root_path,'protein_design','Correct_MJ99_Argentina_PUBLICATION','prok250.dat')
# dat_prok250 = pd.read_csv(prok250_old_KZ_fname)
# dat_prok250.rename(columns={'Temp':'topt'}, inplace=True)
# #########################
#########################
#########################


def label(aa,rr,pp):
    label = '%s: '%aa
    if pp<0.001:
        label+= '$R=%.2f^{***}$ '%rr
    elif pp<0.05:
        label+= '$R=%.2f^{**}$ '%rr
    else:
        label+= '$R=%.2f^{*}$ '%rr
    return label



def plot_palette(loc_dat,fname='aa_palette',vmin=24,vmax=63):
    xbins = 5
    ybins = 4
    # create grid of subplots on a large figure canvas
    # share x&y axes among all plots
    fig, ax = plt.subplots(ybins, xbins, figsize=(7.5,7.5*5.0*0.25), sharex=True, sharey=True)
    # no space between subplots
    fig.subplots_adjust(hspace=0.0, wspace=0.0)
    # make some room for the axes' labels
    l,b,r,t = 0.03,0.3,0.99,0.99
    fig.subplots_adjust(bottom=b, left=l, right=r, top=t)
    w = r-l
    h = t-b
    ##############################
    # lower axes panel for amino acid combinations ...
    w = w/float(xbins)
    height = h = h/float(ybins)
    palette_bottom = b
    bottom = b = 0.04
    right = r
    ax_comb = []
    for i_ax in xrange(xbins):
        ax_comb.append(fig.add_axes([l,b,w,h]))
        if i_ax:
            ax_comb[i_ax].yaxis.set_tick_params(which='both',labelleft='off')
        ax_comb[i_ax].xaxis.set_major_locator( MaxNLocator(nbins = 5,prune='upper') )
        ax_comb[i_ax].yaxis.set_major_locator( MaxNLocator(nbins = 5) )
        ax_comb[i_ax].set_ylim((18,49))
        ax_comb[i_ax].set_xlim((0,115))
        # else:
        #     ax_comb[i_ax].xaxis.set_major_locator( MaxNLocator(nbins = 5,prune='upper') )
        l += w
    ##############################
    # temperature range for the plots:
    t_range = pd.np.asarray([loc_dat['topt'].min(),loc_dat['topt'].max()])
    t_range = (t_range - t_range.mean())*1.1 + t_range.mean()
    #
    ##############################
    # plot combinations ...
    aa_combinations = ['IVYWREL', 'DEKR', 'AGNQSTHY', 'MPCLVWIF', 'ILVM']
    comb_to_plot = [loc_dat[list(aa_comb)].sum(axis=1) for aa_comb in aa_combinations]
    for axx,aa_comb,comb_dat in zip(ax_comb,aa_combinations,comb_to_plot):
        axx.scatter(loc_dat['topt'],comb_dat,edgecolor='none',s=65,c=loc_dat['GC'],vmin=vmin,vmax=vmax,cmap = plt.get_cmap('jet'))
        a,b,r,pval,_ = st.linregress(loc_dat['topt'],comb_dat)
        axx.plot(t_range,a*t_range+b,'-',color='dimgray',lw=2,label=label(aa_comb,r,pval))
        #
        axes_equator = pd.np.asarray(axx.get_ylim()).mean()
        loc = (0.01,0.91) if comb_dat.mean()<axes_equator else (0.01,0.15)
        axx.text(loc[0],loc[1],label(aa_comb,r,pval),fontsize=8.3,fontweight='bold',verticalalignment='top',transform=axx.transAxes)
    ##############################
    for yi in xrange(ybins):
        for xi in xrange(xbins):
            # figuring out corresponding amino acid ...
            aa_num = yi*xbins + xi
            aa = aacids[aa_num]
            # x&y data for plotting in a given axis ...
            scatter = ax[yi,xi].scatter(loc_dat['topt'],loc_dat[aa],edgecolor='none',s=65,c=loc_dat['GC'],vmin=vmin,vmax=vmax,cmap = plt.get_cmap('jet'))
            a,b,r,pval,_ = st.linregress(loc_dat['topt'],loc_dat[aa])
            ax[yi,xi].plot(t_range,a*t_range+b,'-',color='dimgray',lw=2,label=label(aa,r,pval))
            #
            ax[yi,xi].set_ylim((0,14))
            axes_equator = pd.np.asarray(ax[yi,xi].get_ylim()).mean()
            loc = (0.081,0.9) if loc_dat[aa].mean()<axes_equator else (0.081,0.15)
            ax[yi,xi].text(loc[0],loc[1],label(aa,r,pval),fontsize=12,fontweight='bold',verticalalignment='top',transform=ax[yi,xi].transAxes)
            #
    # 
    ax[ybins-2,0].yaxis.set_major_locator( MaxNLocator(nbins = 5,prune='upper') )
    ax[0,xbins-2].xaxis.set_major_locator( MaxNLocator(nbins = 5,prune='upper') )
    # #
    #
    ax[0,xbins-2].set_xlim((0,115))
    #
    # fig.text(0.03,0.5,'composition, %',rotation='vertical',transform=fig.transFigure,fontsize=15,ha='center',va='center')
    fig.text(0.5,0.01,'Temperature, $^{o}C$',transform=fig.transFigure,fontsize=15,ha='center',va='center')
    # #
    left = 0.5
    cax_height = 0.028
    cax_bottom = 0.5*(bottom+height+palette_bottom)-0.5*cax_height
    cax = fig.add_axes([left,cax_bottom,right-left,cax_height])
    fig.text(left-0.2, cax_bottom+0.5*cax_height,'GC content, %',transform=fig.transFigure,fontsize=14,ha='left',va='center')
    cbar = fig.colorbar(scatter,cax=cax,orientation='horizontal')
    ticks = 5*(pd.np.arange(vmin//5,vmax//5)+1)
    ticklabels = map(str,ticks)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(ticklabels)
    # #
    fig.savefig(os.path.join(results_path,'%s.pdf'%fname))




plot_palette(dat)
# plot_palette(dat_bact,fname='aa_palette_bacter',vmin=24,vmax=75)




#################################################
#################################################
#################################################
#################################################
def dataset_quality_plot(x,y,xmin=24,xmax=75,ymin=0,ymax=115,fname='dataset_plot',ylabel='Temperature, $^{o}C$',xlabel='GC content, %'):
    # plt.clf()
    nullfmt = NullFormatter() # no labels
    ###############################
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    #
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    #
    # start with a rectangular Figure
    plt.clf()
    fig = plt.figure(figsize=(7.5,1.0*7.5))
    # add axes ...
    axScatter = fig.add_axes(rect_scatter)
    axHistx = fig.add_axes(rect_histx)
    axHisty = fig.add_axes(rect_histy)
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    #
    # the scatter plot:
    axScatter.scatter(x, y,s=65, edgecolor='none')
    # now determine nice limits by hand:
    num_points = 20
    xbins = pd.np.linspace(xmin,xmax,num=num_points)
    ybins = pd.np.linspace(ymin,ymax,num=num_points)
    #
    axScatter.set_xlim( (xmin, xmax) )
    axScatter.set_ylim( (ymin, ymax) )
    #
    axScatter.set_xlabel(xlabel,fontsize=14)
    axScatter.set_ylabel(ylabel,fontsize=14)
    #
    axHistx.hist(x, bins=xbins,edgecolor='none')
    axHisty.hist(y, bins=ybins,edgecolor='none', orientation='horizontal')
    #
    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )
    #
    axHistx.yaxis.set_major_locator( MaxNLocator(nbins = 5) )
    axHisty.yaxis.set_major_locator( MaxNLocator(nbins = 5) )
    axHisty.xaxis.set_major_locator( MaxNLocator(nbins = 5) )
    #
    # plt.show()
    fig.savefig(os.path.join(results_path,'%s.pdf'%fname))



dataset_quality_plot(dat.GC,dat.topt,xmin=24,xmax=63,ymin=0,ymax=115,fname='dataset_plot')
# dataset_quality_plot(dat_bact.GC,dat_bact.topt,xmin=24,xmax=75,ymin=0,ymax=115,fname='bact_dataset_plot')




# correlations_table = {}
# correlations_table['Prokaryotes'] = []
# correlations_table['Eubacteria'] = []
# correlations_table['Archaea'] = []
# correlations_table['GC50'] = []
# correlations_table['GC30'] = []
# correlations_table['Methanococci'] = []
# for aa in aacids:
#     a,b,r,pval,_ = st.linregress(dat_prok250['topt'],dat_prok250[aa])
#     correlations_table['Prokaryotes'].append(r)
#     a,b,r,pval,_ = st.linregress(dat_bact['topt'],dat_bact[aa])
#     correlations_table['Eubacteria'].append(r)
#     a,b,r,pval,_ = st.linregress(dat[indcs]['topt'],dat[indcs][aa])
#     correlations_table['Archaea'].append(r)
#     a,b,r,pval,_ = st.linregress(dat[indcs_GC50]['topt'],dat[indcs_GC50][aa])
#     correlations_table['GC50'].append(r)
#     a,b,r,pval,_ = st.linregress(dat[indcs_GC30]['topt'],dat[indcs_GC30][aa])
#     correlations_table['GC30'].append(r)
#     a,b,r,pval,_ = st.linregress(dat[indcs_methan]['topt'],dat[indcs_methan][aa])
#     correlations_table['Methanococci'].append(r)

# datasets = ['Prokaryotes', 'Eubacteria', 'Archaea', 'GC50', 'GC30', 'Methanococci']


# correlations_table = pd.DataFrame(correlations_table)
# corr_tab_df = correlations_table.set_index(pd.np.asarray(aacids))
# corr_tab_df = corr_tab_df[datasets]

# plt.clf()
# fig = plt.figure(figsize=(2.05,3.7))
# l,b = 0.15, 0.18
# w,h = 0.92-l, 0.97-b
# ax = fig.add_axes([l,b,w,h])
# # fig, ax = plt.subplots(1,1, figsize=(3.5,3.3))

# # ax.imshow(corr_tab_df.get_values(),interpolation='nearest')
# # image = ax.pcolor(corr_tab_df.get_values()*10.0,vmin=-10,vmax=10,cmap='seismic',norm = mpl.colors.SymLogNorm(linthresh=2.0,linscale=0.0001))
# image = ax.pcolor(corr_tab_df.get_values(),vmin=-1,vmax=1,cmap='seismic')
# ax.set_xticks(pd.np.arange(corr_tab_df.columns.size)+0.5)
# ax.set_xticklabels(corr_tab_df.columns,rotation='vertical')
# ax.set_yticks(pd.np.arange(corr_tab_df.index.size)+0.5)
# ax.set_yticklabels(corr_tab_df.index)
# fig.colorbar(image,ax=ax,orientation='vertical')

# fig.savefig(os.path.join(results_path,'table.pdf'))






































