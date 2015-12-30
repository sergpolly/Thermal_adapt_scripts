import os
import sys
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats as st
import random as rnd
import numpy as np
#################################
# before we proceed to plotting, add the TrOp status calculator for organisms ...
def get_one_trop(all_cds_grouped,idx):
    org_cds = all_cds_grouped.get_group(idx)
    # check if TrOp ...
    # for a given organism(id) all TrOp values must be same
    trop_vals = org_cds['TrOp'].unique()
    assert trop_vals.size == 1
    # then just figure out TrOp value after unpacking ...
    trop, = trop_vals
    if pd.isnull(trop):
        # special return - not enough ribosomal proteins ...
        return 'none'
    if not trop:
        # False, return False 
        return 'false'
    elif trop == True:
        # if it's True just return ...
        return 'true'
    else:
        raise ValueError
#################################
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullFormatter
import scipy.interpolate as interpol
font = {'family' : 'sans-serif',
        #'weight' : 'bold',
        'size'   :9}
mpl.rc('font', **font)
###########################################
aacids = list('CMFILVWYAGTSNQDEHRKP')
aa_combinations = ['IVYWREL', 'DEKR', 'AGNQSTHY', 'MPCLVWIF', 'ILVM']
##########################################
root_path = os.path.expanduser('~')
bact_path = os.path.join(root_path,'GENOMES_BACTER_RELEASE69/genbank')
arch_path = os.path.join(root_path,'GENOMES_ARCH_SEP2015')
# SOME ARCHAEAL DATA ...
arch        = pd.read_csv(os.path.join(arch_path,'summary_organisms_interest.dat'))
arch_nohalo = pd.read_csv(os.path.join(arch_path,'summary_organisms_interest_no_halop.dat'))
###########################################
# SOME BACTERIAL DATA ...
# complete genomes only ...
bact        = pd.read_csv(os.path.join(bact_path,'env_catalog_compgenome.dat'))
# bacter proteomic summary ...
bact_prot = pd.read_csv(os.path.join(bact_path,'proteome_all.dat'))
# arch proteomic summary ...
arch_prot = pd.read_csv(os.path.join(arch_path,'proteome_arch.dat'))
arch_prot[aacids] = arch_prot[aacids]*100.0
#
arch_dat = pd.merge(arch,arch_prot,on='assembly_accession')
arch_nohalo_dat = pd.merge(arch_nohalo,arch_prot,on='assembly_accession')
arch_halo_dat = arch_dat[~arch_dat['assembly_accession'].isin(arch_nohalo['assembly_accession'])]
#
bact_dat = pd.merge(bact,bact_prot,on='GenomicID')
bact_dat[aacids] = bact_dat[aacids]*100.0


calculate_TrOp = True
if calculate_TrOp:
    # we need the following to calculate TrOp status on the fly ...
    ###############################################
    # complete_CDS_CAI_DNA.dat same thing ...
    arch_cai_fname = os.path.join(arch_path,"complete_arch_CDS_CAI_DNA.dat")
    bact_cai_fname = os.path.join(bact_path,"complete_CDS_CAI_DNA.dat")
    #
    arch_cai = pd.read_csv(arch_cai_fname)
    bact_cai = pd.read_csv(bact_cai_fname)
    #
    bact_cai_by_org = bact_cai.groupby('GenomicID')
    arch_cai_by_org = arch_cai.groupby('assembly_accession')
    #
    # calculate TrOp
    bact_dat['TrOp'] = [get_one_trop(bact_cai_by_org,idx) for idx in bact_dat['GenomicID']]
    arch_nohalo_dat['TrOp'] = [get_one_trop(arch_cai_by_org,idx) for idx in arch_nohalo_dat['assembly_accession']]





############################################################
#  PLOTTING FUNCTIONS ...
##############################################################

def aap_get_axes(ranges):
    #
    # limit exapnsion function ...
    def exp_lim(xlims,expansion_coeff = 1.1):
        xmean = 0.5*sum(xlims)
        half_xdelta = 0.5*(xlims[1] - xlims[0])
        new_delta = expansion_coeff*half_xdelta
        return (xmean - new_delta, xmean + new_delta)
    # get ranges before you define axis ...
    # exctracting limits ...
    # and adjusting them a bit ...
    exp_coeffs = (1.1, 1.1, 1.0, 1.1)
    flims, tlims, gclims, fclims = [ exp_lim(lims,coeff) for lims,coeff in zip(ranges,exp_coeffs) ]
    #
    #
    xbins = 5
    ybins = 4
    # create grid of subplots on a large figure canvas
    # share x&y axes among all plots
    fig, ax = plt.subplots(ybins, xbins, figsize=(7.5,xbins*7.5/ybins), sharex=True, sharey=True)
    # no space between subplots
    fig.subplots_adjust(hspace=0.0, wspace=0.0)
    # make some room for the axes' labels
    l,b,r,t = 0.05,0.3,0.98,0.98
    w, h = r-l, t-b
    fig.subplots_adjust(bottom=b, left=l, right=r, top=t)
    #
    # axis lims are predefined by the ranges ...
    # assign just one of them, sharex,sharey will take care of the rest ...
    ax[0,0].set_xlim(tlims)
    ax[0,0].set_ylim(flims)
    #
    ax[-2,0].yaxis.set_major_locator( MaxNLocator(nbins = 5,prune='upper') )
    ax[0,-2].xaxis.set_major_locator( MaxNLocator(nbins = 5,prune='upper') )
    #
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
        # this ylims must be changed further on ...
        ax_comb[i_ax].set_ylim(fclims)
        ax_comb[i_ax].set_xlim(tlims)
        l += w
    #################
    #
    left = 0.5
    cax_height = 0.028
    cax_bottom = 0.5*(bottom+height+palette_bottom)-0.5*cax_height
    cax = fig.add_axes([left,cax_bottom,right-left,cax_height])
    #
    # return axes ...
    return (fig,ax,ax_comb,cax)



def update_ranges(dat,flims=False,tlims=False,gclims=False,fclims=False,temp='OptimumTemperature',gc='GC'):
    # lim of all sort is supposed to be either (min,max) or [min,max] ...
    # fclims is {key1:[min,max], key2:[min,max], ...}
    get_lims = lambda dat_array: (dat_array.min(),dat_array.max())
    #
    def update_lims(new_dat,old_lims):
        # function either return (min,max) from new data, or updates old_lims with the new_dat ...
        if not old_lims:
            return get_lims(new_dat)
        else:
            new_lims = get_lims(new_dat)
            return ( min(new_lims[0],old_lims[0]), max(new_lims[1],old_lims[1]) )
    # we are searching min/max across all 20 amino acids ...
    updated_flims = flims
    for aa in aacids:
        updated_flims = update_lims( dat[aa], updated_flims )
    # yet they share their's organisms GC and Temp ...
    updated_tlims   = update_lims( dat[temp], tlims )
    updated_gclims  = update_lims( dat[gc], gclims )
    #
    # aa combinations update must be done as well ...
    # aa_combinations = ['IVYWREL', 'DEKR', 'AGNQSTHY', 'MPCLVWIF', 'ILVM']
    updated_fclims = fclims
    for combo in aa_combinations:
        updated_fclims = update_lims( dat[list(combo)].sum(axis=1), updated_fclims )
    #
    return (updated_flims, updated_tlims, updated_gclims, updated_fclims)




def fill_palette(dat,axis,coloring_type='solid',color_vlims=False,temp='OptimumTemperature',gc='GC',**kwargs):
    #
    #
    def get_coloring(coloring_type,**kwargs):
        if coloring_type == 'solid':
            if 'color' in kwargs:
                coloring = kwargs['color']
            else:
                coloring = 'red'
        elif coloring_type == 'map':
            coloring = dat[gc]
            # color_vlims must be provided ...
        else:
            raise ValueError("coloring type can either 'solid' or 'map'!")
        # returning ...
        return coloring
    #
    #
    def flabel(aa,rr,pp):
        if bool(rr) and bool(pp):
            label = '%s: '%aa
            if pp<0.001:
                label+= '$R=%.2f^{***}$ '%rr
            elif pp<0.05:
                label+= '$R=%.2f^{**}$ '%rr
            else:
                label+= '$R=%.2f^{*}$ '%rr
            return label
        else:
            return '%s'%aa
    #
    #
    def plot_single_scatter(x,y,ax,flabel,coloring='red',vlims=False,fit=True,alpha=1.0):
        scatter_size = 65
        edgecolor = 'none'
        cmap = plt.get_cmap('jet')
        #
        x_range = np.asarray([x.min(),x.max()])
        # use palette limits or not ...
        kwargs = {'vmin':vlims[0],'vmax':vlims[1]} if vlims else {'norm':True}
        kwargs['alpha'] = alpha
        # plot the scatter ...
        scatter = ax.scatter(x,y,edgecolor=edgecolor,s=scatter_size,c=coloring,cmap=cmap,**kwargs)
        if fit and flabel:
            a,b,r,pval,_ = st.linregress(x,y)
            label = flabel(r,pval)
            ax.plot(x_range,a*x_range+b,'-',color='dimgray',lw=2,label=label)
        elif flabel:
            label = flabel(False,False)
        else:
            label = ""
        #
        axes_equator = np.asarray(ax.get_ylim()).mean()
        loc = (0.02,0.91) if y.mean()<axes_equator else (0.02,0.15)
        ax.text(loc[0],loc[1],label,fontsize=8.3,fontweight='bold',verticalalignment='top',transform=ax.transAxes)
        #
        return scatter
    #
    ###################################################
    # extracting axis ...
    fig,ax,ax_comb,cax = axis
    ybins,xbins = ax.shape
    #
    coloring = get_coloring(coloring_type,**kwargs)
    # vlims will be ignored if coloring is solid ...
    alpha   = kwargs['alpha'] if ('alpha' in kwargs) else 1.0
    fit     = kwargs['fit']   if ('fit'   in kwargs) else True
    label   = kwargs['label'] if ('label' in kwargs) else True
    #
    x = dat[temp]
    for axc,combo in zip(ax_comb,aa_combinations):
        y = dat[list(combo)].sum(axis=1)
        # function: partial argument substitution ...
        label = (lambda rr,pp: flabel(combo,rr,pp)) if label else False
        plot_single_scatter(x,y,axc,flabel=label,coloring=coloring,vlims=color_vlims,fit=fit,alpha=alpha)
    ##############################
    for yi in xrange(ybins):
        for xi in xrange(xbins):
            # figuring out corresponding amino acid ...
            aa_num = yi*xbins + xi
            aa = aacids[aa_num]
            y = dat[aa]
            label = (lambda rr,pp: flabel(aa,rr,pp)) if label else False
            # x&y data for plotting in a given axis ...
            scatter = plot_single_scatter(x,y,ax[yi,xi],flabel=label,coloring=coloring,vlims=color_vlims,fit=fit,alpha=alpha)
    #
    # #
    return scatter


##########################################################
# UNDER CONSTRUCTION ...
##########################################################
def figure_level_caps(axis,scatters,coloring_type,color_vlims=None,scnames=None):
    # extract axis ...
    fig,ax,ax_comb,cax = axis
    # some text for x&y axis ...
    fig.text(0.015,0.5,'composition, %',rotation='vertical',transform=fig.transFigure,fontsize=13,ha='center',va='center')
    fig.text(0.5,0.01,'Temperature, $^{o}C$',transform=fig.transFigure,fontsize=13,ha='center',va='center')
    #
    #
    if (coloring_type == 'map') and (len(scatters)>1):
        raise ValueError("color map supports only 1 set of scatter plots ...")
    elif coloring_type == 'map':
        # extract the only scatter ...
        scatter, = scatters
        # draw that color map (likely GC) ...
        cax.set_visible(True)
        pos = cax.get_position()
        left = pos.x0
        cax_bottom = pos.y0
        cax_height = pos.height
        fig.text(left-0.2, cax_bottom+0.5*cax_height,'GC content, %',transform=fig.transFigure,fontsize=14,ha='left',va='center')
        cbar = fig.colorbar(scatter,cax=cax,orientation='horizontal')
        ticks = 5*np.arange(color_vlims[0]//5,color_vlims[1]//5)+5
        ticklabels = [str(_) for _ in ticks]
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(ticklabels)
    elif coloring_type == 'solid':
        cax.set_visible(False)
        pos = cax.get_position()
        left = pos.x0
        right = left + pos.width
        cax_bottom = pos.y0
        cax_height = pos.height
        fig.legend(scatters, scnames,
                    loc = 'center right',
                    bbox_to_anchor=(right, cax_bottom+0.5*cax_height),
                    bbox_transform=fig.transFigure,
                    scatterpoints=1,
                    markerscale=1.5,
                    ncol = len(scatters))
    #
##########################################################
# UNDER CONSTRUCTION ...
##########################################################


# # that's a working exmaple of closure here ...
# def outer():
#     # defaut definitions ...
#     xlims = [45,45]
#     ylims = [50,50]
#     def inner(x,y):
#         xlims[0],xlims[1] = min(xlims[0],min(x)), max(xlims[1],max(x))
#         ylims[0],ylims[1] = min(ylims[0],min(y)), max(ylims[1],max(y))
#         return (xlims,ylims)
#     return inner
# f = outer()


# drawing functional style ...
def dq_get_axes():
    pass
    # drawing procedures ...
    nullfmt = NullFormatter() # no labels
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    #
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    # start with a rectangular Figure
    plt.clf()
    fig = plt.figure(figsize=(7.5,1.0*7.5))
    # add axes ...
    axScatter   = fig.add_axes(rect_scatter)
    axHistx     = fig.add_axes(rect_histx)
    axHisty     = fig.add_axes(rect_histy)
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    #
    axHistx.yaxis.set_major_locator( MaxNLocator(nbins = 5) )
    axHisty.yaxis.set_major_locator( MaxNLocator(nbins = 5) )
    axHisty.xaxis.set_major_locator( MaxNLocator(nbins = 5) )
    #
    return (fig, axScatter, axHistx, axHisty)



# that's a working exmaple of closure here ...
def dq_stack_data():
    # defaut definitions ...
    xlims = [45,45]
    ylims = [50,50]
    def inner(x,y):
        xlims[0],xlims[1] = min(xlims[0],min(x)), max(xlims[1],max(x))
        ylims[0],ylims[1] = min(ylims[0],min(y)), max(ylims[1],max(y))
        return (xlims,ylims)
    return inner
#


def dq_plot_data(x,y,color,label,xlims,ylims,fig,axScatter,axHistx,axHisty,num_points=30):
    # the scatter plot:
    axScatter.scatter(x,y,s=65, edgecolor='none',facecolor=color,label=label,alpha=0.8)
    # now determine nice limits by hand:
    # num_points = 20
    xbins = np.linspace(*xlims,num=num_points)
    ybins = np.linspace(*ylims,num=num_points)
    #
    axScatter.set_xlim( *xlims )
    axScatter.set_ylim( *ylims )
    #
    axScatter.set_xlabel('GC',fontsize=14)
    axScatter.set_ylabel('Temperature',fontsize=14)
    #
    axHistx.hist(np.asarray(x), bins=xbins,edgecolor='none',facecolor=color, alpha=0.8, log=False, normed=True)
    axHisty.hist(np.asarray(y), bins=ybins,edgecolor='none',facecolor=color, orientation='horizontal',alpha=0.8, log=False, normed=True)
    #
    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )
    axScatter.legend(loc=(0.98,1.03),scatterpoints=1,fontsize=11,handletextpad=0.0,frameon=False)


###########################################
#  ACTUAL PLOTTING USING DEFINED FUNCS ...
############################################


###############################################################
# Supplementary Figure 1.
# constructing and drawing the Dataset Quality plot ...
axes = dq_get_axes()
update_lims = dq_stack_data()
lims = update_lims(arch_halo_dat['GC'],arch_halo_dat['OptimumTemperature'])
lims = update_lims(arch_nohalo_dat['GC'],arch_nohalo_dat['OptimumTemperature'])
lims = update_lims(bact_dat['GC'],bact_dat['OptimumTemperature'])
#
dq_plot_data(arch_halo_dat['GC'],arch_halo_dat['OptimumTemperature'],'red','Halophiles',*(lims+axes))
dq_plot_data(arch_nohalo_dat['GC'],arch_nohalo_dat['OptimumTemperature'],'blue','Archaea',*(lims+axes))
dq_plot_data(bact_dat['GC'],bact_dat['OptimumTemperature'],'green','Bacteria',*(lims+axes))
plt.savefig("SuppFig1.pdf")



###############################################################
# Supplementary Figure 1. (TrOp vs All for bacter)
# constructing and drawing the Dataset Quality plot ...
axes = dq_get_axes()
update_lims = dq_stack_data()
# lims = update_lims(arch_halo_dat['GC'],arch_halo_dat['OptimumTemperature'])
# lims = update_lims(arch_nohalo_dat['GC'],arch_nohalo_dat['OptimumTemperature'])
dat_one = bact_dat[bact_dat['TrOp']=='true']
dat_two = bact_dat[bact_dat['TrOp']!='true']
lims = update_lims(dat_one['GC'],dat_one['OptimumTemperature'])
lims = update_lims(dat_two['GC'],dat_two['OptimumTemperature'])
#
dq_plot_data(dat_one['GC'],dat_one['OptimumTemperature'],'red','CUS Bacteria',*(lims+axes))
dq_plot_data(dat_two['GC'],dat_two['OptimumTemperature'],'blue','Non-CUS Bacteria',*(lims+axes))
plt.savefig("SuppFig1_bact_trop.pdf")



###############################################################
# Supplementary Figure 1. (TrOp vs All for arch)
# constructing and drawing the Dataset Quality plot ...
axes = dq_get_axes()
update_lims = dq_stack_data()
# lims = update_lims(arch_halo_dat['GC'],arch_halo_dat['OptimumTemperature'])
# lims = update_lims(arch_nohalo_dat['GC'],arch_nohalo_dat['OptimumTemperature'])
dat_one = arch_nohalo_dat[arch_nohalo_dat['TrOp']=='true']
dat_two = arch_nohalo_dat[arch_nohalo_dat['TrOp']!='true']
lims = update_lims(dat_one['GC'],dat_one['OptimumTemperature'])
lims = update_lims(dat_two['GC'],dat_two['OptimumTemperature'])
#
dq_plot_data(dat_one['GC'],dat_one['OptimumTemperature'],'red','CUS Archaea',*(lims+axes))
dq_plot_data(dat_two['GC'],dat_two['OptimumTemperature'],'blue','Non-CUS Archaea',*(lims+axes))
plt.savefig("SuppFig1_arch_trop.pdf")






###############################################################
# Supplementary Figure 2.(TrOp vs All for bacter)
dat_one = bact_dat[bact_dat['TrOp']!='true']
dat_two = bact_dat[bact_dat['TrOp']=='true']
#
ranges = update_ranges(dat_one)
# ranges = update_ranges(arch_nohalo_dat,*ranges)
ranges = update_ranges(dat_two,*ranges)
axis = aap_get_axes(ranges)
# only GC lims are needed for "fill_palette" ...
scatter1 = fill_palette(dat_one,axis,coloring_type='solid',color_vlims=ranges[2],temp='OptimumTemperature',gc='GC',color='green',alpha=0.7)
scatter2 = fill_palette(dat_two,axis,coloring_type='solid',color_vlims=ranges[2],color='red',alpha=0.7,fit=False,label=False)
# ...
scnames = list(reversed(['Tr.Op.Bacteria','Non-Tr.Op.Bacteria']))
figure_level_caps(axis,scatters=[scatter1,scatter2],coloring_type='solid',color_vlims=ranges[2],scnames=scnames)
plt.savefig('SuppFig2_bact_trop.pdf')




###############################################################
# Supplementary Figure 2.(TrOp vs All for bacter)
dat_one = arch_nohalo_dat[arch_nohalo_dat['TrOp']=='true']
dat_two = arch_nohalo_dat[arch_nohalo_dat['TrOp']!='true']
#
ranges = update_ranges(dat_one)
# ranges = update_ranges(arch_nohalo_dat,*ranges)
ranges = update_ranges(dat_two,*ranges)
axis = aap_get_axes(ranges)
# only GC lims are needed for "fill_palette" ...
scatter1 = fill_palette(dat_one,axis,coloring_type='solid',color_vlims=ranges[2],temp='OptimumTemperature',gc='GC',color='green',alpha=0.7)
scatter2 = fill_palette(dat_two,axis,coloring_type='solid',color_vlims=ranges[2],color='red',alpha=0.7,fit=False,label=False)
# ...
scnames = list(reversed(['Tr.Op.Archaea','Non-Tr.Op.Archaea']))
figure_level_caps(axis,scatters=[scatter1,scatter2],coloring_type='solid',color_vlims=ranges[2],scnames=scnames)
plt.savefig('SuppFig2_arch_trop.pdf')




##################################################
# Supplementary Figure 2.
ranges = update_ranges(arch_nohalo_dat)
# ranges = update_ranges(arch_nohalo_dat,*ranges)
ranges = update_ranges(arch_halo_dat,*ranges)
axis = aap_get_axes(ranges)
# only GC lims are needed for "fill_palette" ...
scatter1 = fill_palette(arch_nohalo_dat,axis,coloring_type='solid',color_vlims=ranges[2],temp='OptimumTemperature',gc='GC',color='green',alpha=0.7)
scatter2 = fill_palette(arch_halo_dat,axis,coloring_type='solid',color_vlims=ranges[2],color='red',alpha=0.7,fit=False,label=False)
# ...
scnames = ['Archaea','Halophilic Archaea']
figure_level_caps(axis,scatters=[scatter1,scatter2],coloring_type='solid',color_vlims=ranges[2],scnames=scnames)
plt.savefig('SuppFig2_halo.pdf')



###############################################################
# Supplementary Figure 2. BACTERIA
ranges = update_ranges(bact_dat)
axis = aap_get_axes(ranges)
# only GC lims are needed for "fill_palette" ...
scatter = fill_palette(bact_dat,axis,coloring_type='map',color_vlims=ranges[2],temp='OptimumTemperature',gc='GC',color='green',alpha=1.0)
figure_level_caps(axis,scatters=[scatter,],coloring_type='map',color_vlims=ranges[2])
plt.savefig('SuppFig2_bact.pdf')


###############################################################
# Supplementary Figure 2. ARCHAEA
ranges = update_ranges(arch_nohalo_dat)
axis = aap_get_axes(ranges)
# only GC lims are needed for "fill_palette" ...
scatter = fill_palette(arch_nohalo_dat,axis,coloring_type='map',color_vlims=ranges[2],temp='OptimumTemperature',gc='GC',color='green',alpha=1.0)
figure_level_caps(axis,scatters=[scatter,],coloring_type='map',color_vlims=ranges[2])
plt.savefig('SuppFig2_arch.pdf')
















