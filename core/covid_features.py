# This Python script contains functions used for feature finder

import numpy as np
import pandas as pd
import seaborn as sns
import copy
from scipy.stats.stats import pearsonr
import math

import sys,os
from pathlib import Path
home = str(Path.home())
core_dir = home+'/repositories/herg/core/'
conf_dir = core_dir+"/conf"
sys.path.insert(0, core_dir)
sys.path.insert(0, conf_dir)

from descriptor_setup import dnames, dlist
from covid import *
from covid_data import *
from buildmodel import calc_topo_descs

from descriptor_setup import dnames, dlist
from covid import *
from covid_data import *
from buildmodel import calc_topo_descs

from collections import OrderedDict
from rdkit import Chem,DataStructs
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if modifying defaults
import operator
DrawingOptions.bondLineWidth=1.8
IPythonConsole.ipython_useSVG=False

import matplotlib.pyplot as plt


def find_features(in_df, min_range=0.7, max_range=1.0):
    smiles=list(in_df['Structure'])
    mols=[]
    for i in range(len(smiles)):
        mols.append(Chem.MolFromSmiles(smiles[i]))
    descs=calc_topo_descs(mols)

    # number of compounds
    n_mols=descs.shape[0]
    # number of features (default from RDkit)
    n_features=descs.shape[1]
    descs.shape
    descs0=descs

    ## extra backup of original features
    input_descs=descs

    # create a local copy of dnames
    local_dnames = copy.copy(dnames)  
    # array of indices of deleted entries from dnames
    del_indices = []  

    hold = np.copy(input_descs)
    hold.shape

    # Remove zero-variance descriptors
    for i in range((len(local_dnames) - 1), -1, -1):
        if min(input_descs[:, i]) == max(input_descs[:, i]):
            #print(input_descs[:, i]==0)
            if sum(input_descs[:, i])==0:
                hold = np.delete(hold, [i], axis=1)
                del_indices.append(i)
                del local_dnames[i]
    descs = hold
    descs.shape

    # Remove highly correlated descriptors
    # only keep R in certain range
    # also remove nan/inf R
    correl = np.zeros((len(local_dnames), len(local_dnames)))
    hcpairs = []
    hcdescs = {}
    descs_to_del = []

    for i in range(len(local_dnames)):
        for j in range(len(local_dnames)):
            try:
                correl[i, j] = pearsonr(descs[:, i], descs[:, j])[0]
            except:
                correl[i, j] = 0

    for i in range((len(local_dnames) - 1)):
        for j in range(i + 1, len(local_dnames)):
            ################################################################################################
            if (correl[i, j] > min_range) and (correl[i, j] < max_range):
                hcpairs.append((i, j))

    for pair in hcpairs:
        if pair[0] not in hcdescs.keys():
            hcdescs[pair[0]] = 1
        else:
            hcdescs[pair[0]] = hcdescs[pair[0]] + 1
        if pair[1] not in hcdescs.keys():
            hcdescs[pair[1]] = 1
        else:
            hcdescs[pair[1]] = hcdescs[pair[1]] + 1

    sorted_hcdescs = sorted(hcdescs.items(), key=operator.itemgetter(0))
    sorted_hcdescs.reverse()

    while len(sorted_hcdescs) > 0:
        foo = sorted_hcdescs[0][0]
        descs_to_del.append(foo)

        for i in range((len(hcpairs) - 1), -1, -1):
            if foo in hcpairs[i]:
                del hcpairs[i]

        hcdescs = {}
        for pair in hcpairs:
            if pair[0] not in hcdescs.keys():
                hcdescs[pair[0]] = 1
            else:
                hcdescs[pair[0]] = hcdescs[pair[0]] + 1
            if pair[1] not in hcdescs.keys():
                hcdescs[pair[1]] = 1
            else:
                hcdescs[pair[1]] = hcdescs[pair[1]] + 1

        sorted_hcdescs = sorted(hcdescs.items(), key=operator.itemgetter(1))
        sorted_hcdescs.reverse()

    descs_to_del.sort()
    descs_to_del.reverse()

    hold = np.copy(descs)

    del_indices = []
    for i in descs_to_del:
        hold = np.delete(hold, [i], axis=1)
        del_indices.append(i)
        del local_dnames[i]
    descs = hold

    cleaned_descriptors = hold

    # Get complete list of indices of descriptors in use
    indices = []

    for desc in local_dnames:
        i = dnames.index(desc)
        indices.append(i)
    indices.sort()
    
    return indices

def evaluate_features(features_list, col):
    ref = features_list[col]
    overlaps = []
    for features in features_list:
        if features != ref:
            overlaps = set(overlaps) | set(features)
    output = set(ref) - set(overlaps)
    return output

def overlap_features(features_list):
    overlaps = features_list[0]
    for features in features_list:
        overlaps = set(overlaps) & set(features)
    output = set(overlaps)
    return output

def get_descs(in_df):
    smiles=list(in_df['Structure'])
    mols=[]
    for i in range(len(smiles)):
        mols.append(Chem.MolFromSmiles(smiles[i]))
    descs=calc_topo_descs(mols)
    for i in range(descs.shape[0]):
        for j in range(descs.shape[1]):
            if math.isnan(descs[i,j]) or math.isinf(descs[i,j]):
                descs[i,j]=0

    return descs

def loop_features(df, sel_class=['TN','FP']):
    features = []
    for methodsars in lsmethodsars:
        for function in lsfunction:
            for protein in lsprotein:
                df2=get_df_subset(df,
                                  sel_class=sel_class,
                                  sel_sars_method=[methodsars],
                                  sel_er_method=['bla'],
                                  sel_fun=[function],
                                  sel_protein=[protein])
                features.append(find_features(df2))
    return features

def plot_dist(descs):
    # figure parameters
    nrow=40
    ncol=5
    figheight=4*nrow
    figwidth=4.5*ncol

    fig=plt.figure(figsize=(figwidth,figheight))

    for i in range(descs.shape[1]):
        ax = plt.subplot2grid((nrow, ncol), (i//ncol, i%ncol))
        hist, bins = np.histogram(descs[:,i],density=False)
        ## calculate probability instead of using density to aviod sum is not 1 
        hist = hist/sum(hist)

        # manipulate bin
        bin_shift=(bins[1]-bins[0])/2
        nbin_center=len(bins)-1
        xbin=bins+bin_shift
        xbin=xbin[0:nbin_center]

        # smooth x,y datapoints
        xbin_smooth,yhist_smooth=get_smooth(xbin,hist)

        fig.subplots_adjust(hspace = 0.3, wspace=0.4)
        ax.plot(xbin_smooth, yhist_smooth, lw=2, color='blue')
#         ax.hist(list(descs[:,i]), density=False)
        ax.set_xlabel('Distance')
        ax.set_ylabel('Probability')
        ax.set_title(dnames[i])

def compare_dist(descs, viral, antiviral, df_ref, ref):
    # figure parameters
    nrow=40
    ncol=5
    figheight=4*nrow
    figwidth=4.5*ncol

    fig=plt.figure(figsize=(figwidth,figheight))

    for i in range(descs.shape[1]):
        ax = plt.subplot2grid((nrow, ncol), (i//ncol, i%ncol))
        hist_v, bins_v = np.histogram(viral[:,i],density=True)
        hist_av, bins_av = np.histogram(antiviral[:,i],density=True)

        ## calculate probability instead of using density to aviod sum is not 1 
#         hist_v = hist_v/sum(hist_v)
#         hist_av = hist_av/sum(hist_av)


        # manipulate bin
        bin_shift_v=(bins_v[1]-bins_v[0])/2
        nbin_center_v=len(bins_v)-1
        xbin_v=bins_v+bin_shift_v
        xbin_v=xbin_v[0:nbin_center_v]
        
        bin_shift_av=(bins_av[1]-bins_av[0])/2
        nbin_center_av=len(bins_av)-1
        xbin_av=bins_av+bin_shift_av
        xbin_av=xbin_av[0:nbin_center_av]

        # smooth x,y datapoints
        xbin_smooth_v,yhist_smooth_v=get_smooth(xbin_v,hist_v)
        xbin_smooth_av,yhist_smooth_av=get_smooth(xbin_av,hist_av)

        fig.subplots_adjust(hspace = 0.3, wspace=0.4)
        ax.plot(xbin_smooth_v, yhist_smooth_v, lw=2, color='red', label='viral')
        ax.plot(xbin_smooth_av, yhist_smooth_av, lw=2, color='blue', label='antiviral')

#         ax.hist(list(descs[:,i]), density=False)
        ax.set_xlabel('Distance')
        ax.set_ylabel('Probability')
        ax.set_title(dnames[i])
        
        hist, bins = np.histogram(descs[:,i],density=False)
        ## calculate probability instead of using density to aviod sum is not 1 
        hist = hist/sum(hist)

        # manipulate bin
        bin_shift=(bins[1]-bins[0])/2
        nbin_center=len(bins)-1
        xbin=bins+bin_shift
        xbin=xbin[0:nbin_center]

        # smooth x,y datapoints
        xbin_smooth,yhist_smooth=get_smooth(xbin,hist)

        fig.subplots_adjust(hspace = 0.3, wspace=0.4)
        ax.plot(xbin_smooth, yhist_smooth, lw=2, color='yellow', label='all')
#         ax.hist(list(descs[:,i]), density=False)
        
        #ref = [E2_string, tamoxifen_string, raloxifene_string, topotecan_string, testosterone_string]
        dict_color = {E2_string:'gray', tamoxifen_string:'pink', raloxifene_string:'cyan',
                      topotecan_string:'green', testosterone_string:'orange'}
        for j in range(ref.shape[0]):
            ax.axvline(x=ref[j,i], color=dict_color[df_ref['Name'][j]])

def get_descs_matrix(in_df):
    smiles=list(in_df['Structure'])
    mols=[]
    for i in range(len(smiles)):
        mols.append(Chem.MolFromSmiles(smiles[i]))
    descs=calc_topo_descs(mols)
#     for i in range(descs.shape[0]):
#         for j in range(descs.shape[1]):
#             if math.isnan(descs[i,j]) or math.isinf(descs[i,j]):
#                 descs[i,j]=0
    return descs

def do_binning(listin,nbins,bin_range):
    hist, bins = np.histogram(listin,bins=nbins, range=bin_range,density=True)
    bin_shift=(bins[1]-bins[0])/2
    nbin_center=len(bins)-1
    xbin=bins+bin_shift
    xbin=xbin[0:nbin_center]
    return xbin, hist

def get_corr_ind(descs, min_range=0.7, max_range=1.0):
    # number of compounds
    n_mols=descs.shape[0]
    # number of features (default from RDkit)
    n_features=descs.shape[1]
    descs.shape
    descs0=descs

    ## extra backup of original features
    input_descs=descs

    # create a local copy of dnames
    local_dnames = copy.copy(dnames)  
    # array of indices of deleted entries from dnames
    del_indices = []  

    hold = np.copy(input_descs)
    hold.shape
    # Remove zero-variance descriptors
    for i in range((len(local_dnames) - 1), -1, -1):
        if min(input_descs[:, i]) == max(input_descs[:, i]):
            #print(input_descs[:, i]==0)
            if sum(input_descs[:, i])==0:
                hold = np.delete(hold, [i], axis=1)
                del_indices.append(i)
                del local_dnames[i]
    descs = hold
    descs.shape

    # Remove highly correlated descriptors
    # only keep R in certain range
    # also remove nan/inf R
    correl = np.zeros((len(local_dnames), len(local_dnames)))
    hcpairs = []
    hcdescs = {}
    descs_to_del = []

    for i in range(len(local_dnames)):
        for j in range(len(local_dnames)):
            try:
                correl[i, j] = pearsonr(descs[:, i], descs[:, j])[0]
            except:
                correl[i, j] = 0

    for i in range((len(local_dnames) - 1)):
        for j in range(i + 1, len(local_dnames)):
            ################################################################################################
            if (correl[i, j] > min_range) and (correl[i, j] < max_range):
                hcpairs.append((i, j))

    for pair in hcpairs:
        if pair[0] not in hcdescs.keys():
            hcdescs[pair[0]] = 1
        else:
            hcdescs[pair[0]] = hcdescs[pair[0]] + 1
        if pair[1] not in hcdescs.keys():
            hcdescs[pair[1]] = 1
        else:
            hcdescs[pair[1]] = hcdescs[pair[1]] + 1

    sorted_hcdescs = sorted(hcdescs.items(), key=operator.itemgetter(0))
    sorted_hcdescs.reverse()

    while len(sorted_hcdescs) > 0:
        foo = sorted_hcdescs[0][0]
        descs_to_del.append(foo)

        for i in range((len(hcpairs) - 1), -1, -1):
            if foo in hcpairs[i]:
                del hcpairs[i]

        hcdescs = {}
        for pair in hcpairs:
            if pair[0] not in hcdescs.keys():
                hcdescs[pair[0]] = 1
            else:
                hcdescs[pair[0]] = hcdescs[pair[0]] + 1
            if pair[1] not in hcdescs.keys():
                hcdescs[pair[1]] = 1
            else:
                hcdescs[pair[1]] = hcdescs[pair[1]] + 1

        sorted_hcdescs = sorted(hcdescs.items(), key=operator.itemgetter(1))
        sorted_hcdescs.reverse()

    descs_to_del.sort()
    descs_to_del.reverse()

    hold = np.copy(descs)

    del_indices = []
    for i in descs_to_del:
        hold = np.delete(hold, [i], axis=1)
        del_indices.append(i)
        del local_dnames[i]
    descs = hold

    cleaned_descriptors = hold

    # Get complete list of indices of descriptors in use
    indices = []

    for desc in local_dnames:
        i = dnames.index(desc)
        indices.append(i)
    indices.sort()
    
    return indices

def get_correl_matrix(descs):
    # number of compounds
    n_mols=descs.shape[0]
    # number of features (default from RDkit)
    n_features=descs.shape[1]
    descs.shape
    descs0=descs

    ## extra backup of original features
    input_descs=descs

    # create a local copy of dnames
    local_dnames = copy.copy(dnames)  
    # array of indices of deleted entries from dnames
    del_indices = []  

    hold = np.copy(input_descs)
    hold.shape

    descs = hold
    descs.shape

    # Remove highly correlated descriptors
    # only keep R in certain range
    # also remove nan/inf R
    correl = np.zeros((len(local_dnames), len(local_dnames)))
    hcpairs = []
    hcdescs = {}
    descs_to_del = []

    for i in range(len(local_dnames)):
        for j in range(len(local_dnames)):
            try:
                correl[i, j] = pearsonr(descs[:, i], descs[:, j])[0]
            except:
                correl[i, j] = 0
    return correl
                
