"""
Program: covid.py

Written: Kuo Hao Lee
         Andy Guan

All COVID related funcation is here.
"""

import sys
import pathlib

project_path = pathlib.Path('~/repositories/herg')
project_path = project_path.expanduser()
code_dir = project_path / 'common'
sys.path.insert(0, str(code_dir))

from scipy.interpolate import make_interp_spline, BSpline  # smoth function

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import multiprocessing as mp
from collections import OrderedDict
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

## modules for using RDkit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole  # Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, \
    DrawingOptions  # Only needed if modifying defaults

DrawingOptions.bondLineWidth = 1.8
from covid_data import *


## get subset of data based on assay type and class type
def get_subset(df, assay_in, class_in):
    # assay_in = "tox21-erb-bla-p1 vs. CPE"
    # class_in = "TP"
    _df = df[df['Assay'] == assay_in]
    _df2 = _df[_df['class'] == class_in]
    return (_df2)


## test get_subset
# get_subset(stat['Assay'][3],"TP")


## given two dataset
## return total number of molecule and unique number
def get_stat(indf):
    # indf = agonist_cpe_TP

    # total number
    ntot = len(indf['CAS'])
    # unique number
    nuniq = len(indf['CAS'].unique())

    # print(ntot,nuniq)
    return (ntot, nuniq)


## given two dataset
## return number of overlap molecule
def get_overlap(indf1, indf2):
    # indf1 = agonist_cpe_TP
    # indf2 = agonist_cpe_FP
    dfall = [indf1, indf2]
    comp = pd.concat(dfall)
    ntot1, nunique1 = get_stat(indf1)
    ntot2, nunique2 = get_stat(indf2)
    ntot, nunique = get_stat(comp)
    overlap = nunique1 + nunique2 - nunique
    # print(overlap)
    return (overlap)


## given two dataset
## return number of overlap molecule
## method2
def get_overlap2(indf1, indf2):
    # indf1 = agonist_cpe_TP
    # indf2 = agonist_cpe_FP
    overlap = len(list(set(indf1['CAS']) & set(indf2['CAS'])))
    # print(overlap)
    return (overlap)

## given list of dataframes with 1 column (or specify column with key)
## return overlapping elements
## method 3: works independent of having 'CAS' key
def get_overlap3(df_list, key=None):
    if len(df_list) < 1:
        return df_list
    else:
        df = df_list[0]
        i = 1
        if key == None:
            while i < len(df_list):
                df = pd.merge(df, df_list[i], how = 'inner')
                i+=1
        else:
            while i < len(df_list):
                df = pd.merge(df, df_list[i], how = 'inner', on=key, copy=False)
                i+=1    
        return df

## given list of dataframes with 1 column (or specify column with key)
## return nonoverlapping elements in that column
def get_nonoverlap(df_list, key=None):
    if len(df_list) < 1:
        return df_list
    else:
        df = df_list[0]
        i = 1
        if key == None:
            while i < len(df_list):
                comp1 = list(np.setdiff1d(df,df_list[i]))
                comp2 = list(np.setdiff1d(df_list[i],df))
                df = pd.DataFrame(comp1+comp2)
                i+=1
        else:
            while i < len(df_list):
                comp1 = list(np.setdiff1d(df[key],df_list[i][key]))
                comp2 = list(np.setdiff1d(df_list[i][key],df[key]))
                df = pd.DataFrame()
                df[key] = comp1 + comp2
                i+=1    
        return df
    
## given list of conditions and dataframe
## return dataframe with rows containing conditional elements in column @ key
def get_subsetfromconditions(indf, inlist, key):
    df = pd.DataFrame()
    for i in range(len(inlist)):
        df2 = indf.loc[indf[key]==inlist[i]]
        df = pd.concat([df,df2])
    return df

def get_fp(smile1, smile2, fp):
    mol1 = Chem.MolFromSmiles(smile1)
    mol2 = Chem.MolFromSmiles(smile2)
    if fp == fp_morgan:
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, 4096)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, 4096)
    else:
        fp1 = Chem.RDKFingerprint(mol1)
        fp2 = Chem.RDKFingerprint(mol2)
    return fp1, fp2


###==========================================================================
# given two smile string
# return the similarity
# https://rdkit.readthedocs.io/en/latest/GettingStartedInPython.html#fingerprinting-and-molecular-similarity
def get_similarity(smile1, smile2, fp, sim):
    fp1, fp2 = get_fp(smile1, smile2, fp)
    if sim == sim_dice:
        similarity = DataStructs.DiceSimilarity(fp1, fp2)
    else:
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return round(similarity, 3)

# plot individual scatter plot
def plot2d(df2, proteins, ref1, ref2,fp, sim, annotate=False, n_min=2, n_max=20):
    fig=plt.figure(figsize=(14, 3.5 * len(proteins)))
    #fig = plt.figure()
    for iprotein in range(len(proteins)):
        if proteins[iprotein] is None:
            print(ref1 + '-x vs ' + ref2 + '-y')
            protein_string = 'all ER'
        else:
            separator = ' + '
            protein_string = separator.join(proteins[iprotein])
        key1 = simstring + ref1
        key2 = simstring + ref2
        for imethodsars in range(len(lsmethodsars)):
            for ifunction in range(len(lsfunction)):
                ax = plt.subplot2grid((len(proteins), 4), (iprotein, imethodsars * 2 + ifunction))
                for iclass in range(len(lsclass)):
                    df3 = get_df_subset(df2, sel_class=[lsclass[iclass]],
                                        sel_sars_method=[lsmethodsars[imethodsars]],
                                        sel_er_method=['bla'], sel_fun=[lsfunction[ifunction]],
                                        sel_protein=proteins[iprotein])
                    ax.scatter(df3[key1], df3[key2], color=color[iclass], marker=marker[iclass],
                               s=40, label=lsclass[iclass])
                    if iprotein == 0:
                        ax.set_title(lsmethodsars[imethodsars] + "_" + lsfunction[ifunction])
                    if imethodsars == 0 and ifunction == 0:
                        ax.set_ylabel(protein_string)
                    if ifunction > 0 or imethodsars > 0:
                        ax.set_ylabel("")
                        ax.set_yticklabels("")
                    if iprotein < len(proteins) - 1:
                        for tk in ax.get_xticklabels():
                            tk.set_visible(False)
                if (annotate == True) and (n_min>1):
                    df4 = get_df_subset(df2, sel_class=lsclass,
                                        sel_sars_method=[lsmethodsars[imethodsars]],
                                        sel_er_method=['bla'], sel_fun=[lsfunction[ifunction]],
                                        sel_protein=proteins[iprotein])
                    plot_annotations(ax, df4, key1, key2, n_min, n_max)
                ax.set_ylim(-0.02, 1.02)
                ax.set_xlim(-0.06, 1.06)
    plt.tight_layout(pad=0, w_pad=0, h_pad=0)
    #fig = plt.figure()
    fig.text(0.5, -0.02, key1, ha='center')
    fig.text(-0.02, 0.5, key2, va='center', rotation='vertical')
    plt.legend()
    fig.savefig(fp+'-'+sim+'_'+ref1+'-'+ref2+'.png', dpi=300)

def plot2d_nonoverlaps(df2, proteins, ref1, ref2,fp, sim, key, annotate=False, n_min=2, n_max=20):
    fig=plt.figure(figsize=(14, 3.5 * len(proteins)))
    #fig = plt.figure()
    for iprotein in range(len(proteins)):
        if proteins[iprotein] is None:
            print(ref1 + '-x vs ' + ref2 + '-y')
            protein_string = 'all ER'
        else:
            separator = ' + '
            protein_string = separator.join(proteins[iprotein])
        key1 = simstring + ref1
        key2 = simstring + ref2
        for imethodsars in range(len(lsmethodsars)):
            for ifunction in range(len(lsfunction)):
                ax = plt.subplot2grid((len(proteins), 4), (iprotein, imethodsars * 2 + ifunction))
                dflist = []
                for iclass in range(len(lsclass)):
                    dflist.append(get_df_subset(df2, sel_class=[lsclass[iclass]],
                                        sel_sars_method=[lsmethodsars[imethodsars]],
                                        sel_er_method=['bla'], sel_fun=[lsfunction[ifunction]],
                                        sel_protein=proteins[iprotein]))
                uniques = get_nonoverlap(dflist, key)
                uniques = uniques[key].unique()
                df_sel = get_subsetfromconditions(df2, uniques, key)
                for iclass in range(len(lsclass)):
                    df3 = get_df_subset(df_sel, sel_class=[lsclass[iclass]],
                                        sel_sars_method=[lsmethodsars[imethodsars]],
                                        sel_er_method=['bla'], sel_fun=[lsfunction[ifunction]],
                                        sel_protein=proteins[iprotein])
                    ax.scatter(df3[key1], df3[key2], color=color[iclass], marker=marker[iclass],
                               s=40, label=lsclass[iclass])
                    if iprotein == 0:
                        ax.set_title(lsmethodsars[imethodsars] + "_" + lsfunction[ifunction])
                    if imethodsars == 0 and ifunction == 0:
                        ax.set_ylabel(protein_string)
                    if ifunction > 0 or imethodsars > 0:
                        ax.set_ylabel("")
                        ax.set_yticklabels("")
                    if iprotein < len(proteins) - 1:
                        for tk in ax.get_xticklabels():
                            tk.set_visible(False)
                if (annotate == True) and (n_min>1):
                    df4 = get_df_subset(df_sel, sel_class=lsclass,
                                        sel_sars_method=[lsmethodsars[imethodsars]],
                                        sel_er_method=['bla'], sel_fun=[lsfunction[ifunction]],
                                        sel_protein=proteins[iprotein])
                    plot_annotations(ax, df4, key1, key2, n_min, n_max)
                ax.set_ylim(-0.02, 1.02)
                ax.set_xlim(-0.06, 1.06)
    plt.tight_layout(pad=0, w_pad=0, h_pad=0)
    #fig = plt.figure()
    fig.text(0.5, -0.02, key1, ha='center')
    fig.text(-0.02, 0.5, key2, va='center', rotation='vertical')
    plt.legend()
    fig.savefig(fp+'-'+sim+'_'+ref1+'-'+ref2+'.png', dpi=300)

## use k-means clustering to plot outlying data points
def plot_annotations(ax, indf, key1, key2, n_min, n_max):
    df = pd.DataFrame()
    df[key1] = indf[key1]
    df[key2] = indf[key2]
    k_list = []
    for i in range(n_min,n_max):
        kmeans = KMeans(init='random', n_clusters=i, n_init=10, max_iter=300, random_state=42)
        k_list.append(silhouette_score(df,kmeans.fit_predict(df)))
        # TODO: show a plot of number of clusters (i) vs silhouette_score
        # RMSD
    n_clusters = k_list.index(max(k_list))+n_min
    kmeans = KMeans(init='random', n_clusters=n_clusters, n_init=10, max_iter=300, random_state=42)
    # TODO: do this two or more times to see if the results are similar?
    kmeans.fit(df)
    df["cluster"] = kmeans.labels_
    for icluster in range(n_clusters):
        rows = df[df["cluster"] == icluster].index
        #if cluster size is not too big, add annotation of chemical compound name
        if len(rows) <= (3*n_max/n_clusters):
            for irows in range(len(rows)):
                ind = rows[irows]
                label = indf["Name"][ind]
                ax.annotate(s=label,xy = (df[key1][ind],df[key2][ind]))
    
# plot scatter plot according to different protein set
def plot(df2, ref1, ref2, smiles, fp, sim, annotate=False, n_min=2, n_max=20, nonoverlaps=False):
    print('similarity between ' + ref1 + ' and ' + ref2 + ':' + str(
        get_similarity(smiles[ref1], smiles[ref2], fp, sim)))
    proteins = [None, ['err'], ['pgc-err'], ['er'], ['erb']]
    if nonoverlaps == True:
        plot2d_nonoverlaps(df2, proteins, ref1, ref2, fp, sim, "Name", annotate, n_min, n_max)
    else:
        plot2d(df2, proteins, ref1, ref2, fp, sim, annotate, n_min, n_max)
        

def get_sim_to_ref(df, mols_strings, smiles, fp, sim):
    pool = mp.Pool(mp.cpu_count())
    for x in mols_strings:
        if x not in smiles:
            df_sel = df[df['Name'] == x]
            smiles[x] = list(df_sel['Structure'])[0]
        # update dataframe
        key = simstring + x
        if key not in df.columns:
            df[key] = pool.starmap(get_similarity, [(df['Structure'][i], smiles[x], fp, sim) for i
                                                    in range(len(df['Structure']))])
    return df

###==========================================================================

## given reference 2 smile 
## return similiary in sml2ref2 column
def get_update_df(dfin, ref2_smi):
    ## get similarity using ref2
    col_sml2ref2 = []

    for i in range(0, len(dfin['Assay'])):
        # print(i)
        col_sml2ref2.append(get_similarity(dfin['Structure'][i], ref2_smi))

    dfin['sml2ref2'] = col_sml2ref2

    return (dfin)


## get subset of dataframe based on given conditions
def get_df_subset(dfin,
                  sel_class=None, sel_sars_method=None,
                  sel_er_method=None, sel_fun=None, sel_protein=None):
    if type(sel_class) == str:
        sel_class = list(sel_class)
    if type(sel_sars_method) == str:
        sel_sars_method = list(sel_sars_method)
    if type(sel_er_method) == str:
        sel_er_method = list(sel_er_method)
    if type(sel_fun) == str:
        sel_fun = list(sel_fun)
    if type(sel_sars_method) == str:
        sel_protein = list(sel_protein)

    dfsel = dfin

    if str(sel_class) != 'None':
        dfsel = dfsel[dfsel['class'].isin(sel_class)]

    if str(sel_sars_method) != 'None':
        dfsel = dfsel[dfsel['Method_SARS'].isin(sel_sars_method)]

    if str(sel_er_method) != 'None':
        dfsel = dfsel[dfsel['Method_ER'].isin(sel_er_method)]

    if str(sel_fun) != 'None':
        dfsel = dfsel[dfsel['Function'].isin(sel_fun)]

    if str(sel_protein) != 'None':
        dfsel = dfsel[dfsel['Protein'].isin(sel_protein)]

    return (dfsel)


###==========================================================================


###==========================================================================
## make smooth curve
def get_smooth(x, y):
    # 300 represents number of points to make between T.min and T.max
    xnew = np.linspace(x.min(), x.max(), 300)
    # type: BSpline
    spl = make_interp_spline(x, y, k=2)
    power_smooth = spl(xnew)
    return (xnew, power_smooth)


## generate hist
def get_hist_f(datain, ctrl_nbin=20):
    # ctrl_nbin=20
    # datain=list(dfsel['sml2Estradiol'])
    # histogram the raw data
    hist_f, bins = np.histogram(datain, bins=ctrl_nbin, range=[0, 1], density=False)

    # calculate probability instead of using density to aviod sum is not 1 
    hist_f = hist_f / sum(hist_f)

    # manipulate bin
    bin_shift = (bins[1] - bins[0]) / 2
    nbin_center = len(bins) - 1
    xbin = bins + bin_shift
    xbin = xbin[0:nbin_center]

    x1, y1 = get_smooth(xbin, hist_f)
    return (x1, y1)


###==========================================================================


## the following functions (get_method_sars, get_functionm, get_protein,
## get_method_er and get_type) wroten by Andy Guan
## the purpose is to add labels to dataframe
###==========================================================================
def get_method_sars(assay_in=None):
    method_sars = None

    if 'vs. CPE' in assay_in:
        method_sars = 'CPE'

    if 'vs. PP' in assay_in:
        method_sars = 'PP'

    return method_sars


def get_function(assay_in=None):
    function = 'agonist'

    if 'antagonist' in assay_in:
        function = 'antagonist'

    return function


def get_protein(assay_in=None):
    protein = None

    if 'er' in assay_in:
        protein = 'er'

    if 'erb' in assay_in:
        protein = 'erb'

    if 'err' in assay_in:
        protein = 'err'

        if 'pgc' in assay_in:
            protein = 'pgc-err'

    return protein


def get_method_er(assay_in=None):
    method_er = None

    if 'bla' in assay_in:
        method_er = 'bla'

    if 'luc-bg1-4e2' in assay_in:
        # method_er = 'luc-bg1-4e2'
        method_er = 'luc'
    else:
        method_er = 'bla'  # if no luc label, then that's bla

    return method_er


def get_type(assay_in=None):
    typ = None

    if 'p1' in assay_in:
        typ = 'p1'

    if 'p2' in assay_in:
        typ = 'p2'

    return typ


## Add columns to dataframe
def add_info_columns(df):
    df_col_method_sars = []
    df_col_function = []
    df_col_protein = []
    df_col_method_er = []
    df_col_type = []
    if 'Assay' in df:
        for i in range(len(df)):
            df_col_method_sars.append(get_method_sars(df['Assay'][i]))
            df_col_function.append(get_function(df['Assay'][i]))
            df_col_protein.append(get_protein(df['Assay'][i]))
            df_col_method_er.append(get_method_er(df['Assay'][i]))
            df_col_type.append(get_type(df['Assay'][i]))
        df['Method_SARS'] = df_col_method_sars
        df['Function'] = df_col_function
        df['Protein'] = df_col_protein
        df['Method_ER'] = df_col_method_er
        df['Type'] = df_col_type
    return df
###==========================================================================
