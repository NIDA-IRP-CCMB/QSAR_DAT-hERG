## define enviroment
import sys, os
from pathlib import Path
home = str(Path.home())
base_dir = home+'/repositories/herg/hERGvDAT/'
core_dir = base_dir+'/core'
conf_dir = core_dir+'/conf'
sys.path.insert(0, conf_dir)
sys.path.insert(0, core_dir)

import pandas as pd
import numpy as np
import filter_config as config
import math

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import SaltRemover
from rdkit import DataStructs
from molvs import Standardizer



import warnings
warnings.filterwarnings('ignore')

def read_data(file_name, Verbose=False):
    '''
    labels = {0:'pref_name', 1:'organism', 2:'assay_id', 3:'assay_type', 
                   4:'relationship_type', 5:'relationship_desc', 6:'confidence_score',
                   7'curated_by', 8:'description', 9:'activity_id', 10:'relation',
                   11:'value', 12:'units', 13:'type', 14:'standard_relation', 15:'standard_value',
                   16:'standard_units', 17:'standard_flag', 18:'standard_type', 19:'pchembl_value',
                   20:'activity_comment', 21:'data_validity_comment', 22:'potential_duplicate',
                   23:'text_value', 24:'standard_text_value', 25:'molregno', 26:'chembl_id', 
                   27:'canonical_smiles', 28:'pref_name', 29:'parent_molregno', 30:'active_molregno', 
                   31:'doc_id', 32:'pubmed_id', 33:'doi', 34:'journal', 35:'year', 36:'volume', 
                   37:'first_page', 38:'src_short_name'}
    '''
    labels = ['pref_name_target', 'organism', 'assay_id', 'assay_type', \
               'relationship_type', 'relationship_desc', 'confidence_score', \
               'curated_by', 'description', 'activity_id', 'relation', \
               'value', 'units', 'type', 'standard_relation', 'standard_value', \
               'standard_units', 'standard_flag', 'standard_type', 'pchembl_value', \
               'activity_comment', 'data_validity_comment', 'potential_duplicate', \
               'text_value', 'standard_text_value', 'molregno', 'chembl_id', \
               'canonical_smiles', 'pref_name', 'parent_molregno', 'active_molregno', \
               'doc_id', 'pubmed_id', 'doi', 'journal', 'year', 'volume', 'first_page', \
               'src_short_name']
    raw_data = pd.read_csv(file_name, names=labels, header=None, sep="\t")

    if Verbose:
        print('Number of compounds at starting: ', len(raw_data))

    return raw_data


def filter_confidence(in_lines, broad=False, kickout=False, Verbose=False):

    # Remove compounds with a confidence score less than 9
    # this should be run early in the filter pipeline

    if broad:
        #one line version below, takes 0.5 seconds (longer)
        #in_lines = in_lines[in_lines['confidence_score'].apply(lambda x: np.any(np.in1d(x, [8,9])))]
        in_lines1=in_lines[in_lines['confidence_score']== 8]
        in_lines2=in_lines[in_lines['confidence_score']== 9]
        in_lines = pd.concat([in_lines1,in_lines2])
    else:
        in_lines=in_lines[in_lines['confidence_score']==9]

    if Verbose:
        print('Number of compounds after confidence score filter: ', len(in_lines))

    return in_lines.reset_index(drop=True)


def filter_assay_type(in_lines, kickout=False, Verbose=False):

    # Remove entries that are not binding or functional studies
    # if this filter is used, it should be done early in the pipeline
    
    #one line version below, takes longer
    #in_lines = in_lines[in_lines['confidence_score'].apply(lambda x: np.any(s in x for s in ['B','F']))]
    
    in_lines1 = in_lines[in_lines['assay_type'] == 'B']
    in_lines2 = in_lines[in_lines['assay_type'] == 'F']
    in_lines = pd.concat([in_lines1,in_lines2])

    if Verbose:
        print('Number of compounds after assay type filter: ', len(in_lines))

    return in_lines.reset_index(drop=True)


def filter_affinity(in_lines, keepKi=True, keepIC50=True, kickout=False, Verbose=False):

    # Remove entries that are not Ki or IC50 values
    # if this filter is used, it should be placed early in the filter pipeline
    
    in_lines1 = pd.DataFrame()
    in_lines2 = pd.DataFrame()
    
    if keepKi:
        in_lines1 = in_lines[in_lines['standard_type'] == 'Ki']
    if keepIC50:
        in_lines2 = in_lines[in_lines['standard_type'] == 'IC50']

    in_lines = pd.concat([in_lines1,in_lines2])

    if Verbose:
        print('Number of compounds after Ki / IC50 filter: ', len(in_lines))

    return in_lines.reset_index(drop=True)


def filter_units(in_lines, kickout=False, Verbose=False):

    # Remove compounds without a standard unit of nM
    # this filter should be placed after the filter_affinity filter is used to limit the number of false positives
    
    in_lines = in_lines[in_lines['standard_units'] == 'nM']

    if Verbose:
        print('Number of compounds after standard units filter: ', len(in_lines))

    return in_lines.reset_index(drop=True)

def filter_exact(in_lines, include_ceilings=False, include_drugmatrix=False, kickout=False, Verbose=False):

    # Process compounds with inexact relationships
    # this should probably be run early in the filtering pipeline, but after values that are not Ki or IC50 are removed
    
    in_lines_sel = in_lines[in_lines['pchembl_value'].isnull()]
    in_lines1 = in_lines[in_lines['standard_relation'] == '=']
    in_lines2 = pd.DataFrame()
    in_lines3 = pd.DataFrame()
    
    #the code in this if statement works, but results in a red message
    if include_ceilings:
        in_lines_sel[in_lines_sel['standard_relation'] == '>']['pchembl_value'] = 1.0
        in_lines_sel[in_lines_sel['standard_relation'] == '>=']['pchembl_value'] = 1.0
        in_lines_sel[in_lines_sel['standard_relation'] == '<']['pchembl_value'] = 11.0
        in_lines_sel[in_lines_sel['standard_relation'] == '<=']['pchembl_value'] = 11.0
        in_lines2 = in_lines_sel[in_lines_sel['standard_relation'] != '=']
        in_lines2 = in_lines2.dropna(thresh=1)
    
    if include_drugmatrix:
        in_lines_sel2 = in_lines_sel[in_lines_sel['src_short_name'] == 'DRUGMATRIX']
        in_lines_sel2[in_lines_sel2.activity_comment.str.startswith('Not Active')]['chembl_value'] = 1.0
        in_lines3 = in_lines_sel2
    
    in_lines = pd.concat([in_lines1,in_lines2,in_lines3])

    if Verbose:
        print('Number of compounds after activity relationship type fixes: ', len(in_lines))

    return in_lines.reset_index(drop=True)

def filter_assaydefinition(in_lines, target, key, kickout=False,
                                Verbose=False):
    # filter for displacement assay data
    filterfile = conf_dir+'/assaydefinition_'+target+'_'+key+'.txt'
    df = pd.read_table(filterfile, names=['keys', 'text'])
    selection = list(df['text'])

    in_lines_in = in_lines[
        in_lines['description'].apply(lambda x: any([s in x for s in selection]))]
    in_lines_out = in_lines[
        in_lines['description'].apply(lambda x: all([s not in x for s in selection]))]
    
    if Verbose:
        print('Number of compounds in ' + key, len(in_lines_in))
        in_lines_in[['description', 'pubmed_id', 'doi']].to_csv("hERG_data_" + key + ".dat",
                                                                sep='\t', index=False)
        in_lines_in.to_csv("hERG_data_" + key + ".tsv", sep='\t', index=False)

        print('Number of compounds out', len(in_lines_out))

    return in_lines_in.reset_index(drop=True), in_lines_out.reset_index(drop=True)


def filter_secondary_test_set(in_lines, kickout=False, Verbose=False):

    # Remove compounds present in secondary test set
    # this filter only applies to andy_hERG models

    # maybe we should break this line into 2 to make it easier to read
    in_lines = in_lines.drop(in_lines[in_lines['chembl_id'].apply(lambda x: any([s == str(x) for s in config.test_set_2_compounds]))].index)
    
    if Verbose:
        print('Number of compounds after removing testset 2 compounds: ', len(in_lines))

    return in_lines.reset_index(drop=True)

def add_doc_cmpd_count(in_lines, kickout=False, Verbose=False):

    # add field to end of data to represent how many compounds
    # are in a given publication (this makes field 39)
    #
    # This should be the first filter run, since it adds a field that other filters might depend upon

    count = []
    
    #reset index before iterating/looping
    in_lines = in_lines.reset_index()
    mol_per_doc_counts = {}
    
    for i in range(len(in_lines)):
        if in_lines['doc_id'][i] not in mol_per_doc_counts.keys():
            mol_per_doc_counts[in_lines['doc_id'][i]] = 1
        else:
            mol_per_doc_counts[in_lines['doc_id'][i]] = mol_per_doc_counts[in_lines['doc_id'][i]] + 1
    
    for i in range(len(in_lines)):
        count.append(mol_per_doc_counts[in_lines['doc_id'][i]])
    
    in_lines['doc_cmpd_count'] = count

    if Verbose:
        print('Number of compounds at adding mols_per_doc count: ', len(in_lines))

    return in_lines.reset_index(drop=True)

def filter_small_sets(in_lines, threshold=5, kickout=False, Verbose=False):

    # Remove compounds that come from sets of less than threshold compounds
    # this filter is generally only applicable for regression models

    # Add number of molecules per document to data first, so that the filter can be done
    in_lines = add_doc_cmpd_count(in_lines,Verbose=False)
    
    #I think this line causes some warning message
    in_lines = in_lines[in_lines['doc_cmpd_count'] > threshold]

    if Verbose:
        print('Number of compounds after data set size filter: ', len(in_lines))

    return in_lines.reset_index(drop=True)

#should compare dataframes later to check if salts were removed properly
def filter_salts(in_lines, kickout=False, Verbose=False):

    # standardize structures and remove salts
    #
    # This should be called before any other filters having to do with molecular structures as it
    # affects both the molecular structure and the molecular weight of many compounds that come out of ChEMBL

    s = Standardizer()
    #salt_file = code_dir / 'Salts.txt'
    salt_file = conf_dir+'/Salts.txt'
    remover = SaltRemover.SaltRemover(defnFilename=salt_file)
    
    for i in range(len(in_lines)):
        mol_in = Chem.MolFromSmiles(in_lines['canonical_smiles'][i])
        mol_out = s.standardize(mol_in)
        smiles_out = Chem.MolToSmiles(remover(mol_out), isomericSmiles=False)
        if '.' in smiles_out:
            in_lines = in_lines.drop(i)
            if kickout:
                kickouts.writerow(in_lines[i])
        else:
            in_lines.loc[i, 'canonical_smiles'] = smiles_out
#             in_lines['canonical_smiles'].replace(i,smiles_out)
#             ## I believe you should just use replace
# The replace function replaces values equal to i with smiles_out
# so I do not think we want to use replace

    if Verbose:
        print('Number of compounds after desalting pass: ', len(in_lines))

    return in_lines.reset_index(drop=True)


def filter_elements(in_lines, kickout=False, Verbose=False):

    # remove entries with oddball elements
    # this needs to be run after the desalting and molecular standardization step
    
    element_filter = Chem.MolFromSmarts('[!C&!c&!N&!n&!O&!o&!S&!s&!P&!p&!F&!Cl&!Br&!I]')
    for i in range(len(in_lines)):
        curr_mol = Chem.MolFromSmiles(in_lines['canonical_smiles'][i])
        if not curr_mol.HasSubstructMatch(element_filter):
            continue
        else:
            in_lines = in_lines.drop(i)
            if kickout:
                kickouts.writerow(in_lines[i])
                
    if Verbose:
        print('Number of compounds after oddball element filter: ', len(in_lines))

    return in_lines.reset_index(drop=True)

def filter_size(in_lines, maxweight=650, kickout=False, Verbose=False):

    # remove compounds with a MW that is greater than the maximum
    # this needs to be run after the structure standardization and desalting step

    for i in range(len(in_lines)):
        molweight = Chem.CalcExactMolWt(Chem.MolFromSmiles(in_lines['canonical_smiles'][i]))
        if molweight >= maxweight:
            in_lines = in_lines.drop(i)

    if Verbose:
        print('Number of compounds after molecular weight filter: ', len(in_lines))

    return in_lines.reset_index(drop=True)
'''
14:'standard_relation', 15:'standard_value',
                   16:'standard_units', 17:'standard_flag', 18:'standard_type', 19:'pchembl_value'
'''
#in progress
def filter_pchembl_values(in_lines, replace=False, kickout=False, Verbose=False):

    # remove compounds that don't have a pChEMBL value associated with them unless we want to
    # calculate it manually.  This should probably be run after filter_exact so that floor and
    # ceiling values get added before they would be eliminated.
    
    nans = in_lines[in_lines['pchembl_value'].isnull()]
    in_lines = in_lines.drop(list(nans.index))
    
    if replace:
        drop_14 = list(nans[nans['standard_relation'] != '='].index)
        drop_15 = list(nans[nans['standard_value'].isnull()].index)
        drops = nans[nans['index'].apply(lambda x: any(s in x for s in 
                                                       list(drop_14 + list(set(drop_15) - set(drop_14)))))]
        nans = nans.drop(list(drops.index))
        nans = nans.reset_index()
        for i in range(len(nans)):
            if nans['standard_units'][i] in ['M','mM','uM','nM','pM','fM']:
                nans.loc[i, 'pchembl_value'] = calc_pscale(nans['standard_value'][i],
                                                           nans['standard_units'][i])
            else:
                drops = pd.concat([drops, nans.iloc[list(i)]])
        nans = nans.drop(list(drops.index), errors = 'ignore')
        in_lines = pd.concat([in_lines,nans])
    else:
        drops = nans
    
    if kickout:
        drops = drops.reset_index(drop=True)
        for i in range(len(drops)):
            kickouts.writerow(drops[i])

    if Verbose:
        print('Number of compounds after pChEMBL value filter: ', len(in_lines))

    return in_lines.reset_index(drop=True)

def filter_weirdos(in_lines, kickout=False, Verbose=False):

    # filter out some odd edge cases we occasionall see in  a full data dump
    # this may not be necessesary to use routinely.

    drop_1 = in_lines[in_lines['standard_value'].isnull()]
    drop_2 = in_lines[in_lines['standard_units'] != 'nM']
    drops = in_lines.iloc[list(drop_1.index + list(set(drop_2.index) - set(drop_1.index)))]
    in_lines = in_lines.drop(drops.index)
    
    if kickout:
        drops = drops.reset_index(drop=True)
        for i in range(len(drops)):
            kickouts.writerow(drops[i])


    if Verbose:
        print('Number of compounds after edge case filter: ', len(in_lines))

    return in_lines.reset_index(drop=True)
'''
0:'pref_name', 1:'organism', 15:'standard_value',
                   16:'standard_units', 17:'standard_flag', 18:'standard_type', 19:'pchembl_value'
'''
def deduplicate_mols(in_lines, limit_range=False, kickout=False, Verbose=False):

    # verify that we only have one instance of a molecule in the final set
    # this needs to be the last filter run

    working_copy = pd.DataFrame()
    fingerprints = []
    mol_dup_pairs = []
    mol_duplicates = []
    holding = []
    seen = []

    # take fingerprints of all remaining molecules

    for i in range(len(in_lines)):
        mol_in = Chem.MolFromSmiles(in_lines['canonical_smiles'][i])
        fp = Chem.GetMorganFingerprintAsBitVect(mol_in, 2, nBits=2048) #set nbits equal to what you will use in model?
        fingerprints.append(fp)

    # pairwise comparison of all fingerprints (matching pairs go into mol_dup_pairs)

    for i in range(0, ((len(fingerprints)) - 1), 1):
        for j in range((i + 1), (len(fingerprints))):
            if DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j]) > 0.999:
                mol_dup_pairs.append((i, j))

    for dup_pair in mol_dup_pairs:
        if dup_pair[0] not in seen and dup_pair[1] not in seen and dup_pair[0] not in holding:
            holding.append(dup_pair[0])
            holding.sort()
            seen.append(dup_pair[0])
            seen.append(dup_pair[1])
            seen.sort()
        if dup_pair[1] not in seen:
            seen.append(dup_pair[1])
            seen.sort()

    for item in holding:
        mol_duplicates.append([item])

    for dup_pair in mol_dup_pairs:
        for i in range(len(mol_duplicates)):
            if dup_pair[0] in mol_duplicates[i] and dup_pair[1] not in mol_duplicates[i]:
                mol_duplicates[i].append(dup_pair[1])

    # select lowest affinity value out of all values per compound (except 1)

    for dup_group in mol_duplicates:

        all_acts = []
        sub_lines = []

        # make certain that we have a pchembl value for each entry and make a list of them
        
        for line_no in dup_group:
            sub_lines.append(in_lines.iloc[[line_no]])
            if math.isnan(in_lines['pchembl_value'][line_no]):
                all_acts.append(calc_pscale(in_lines['standard_value'][line_no],
                                            in_lines['standard_units'][line_no]))
            else:
                all_acts.append(in_lines['pchembl_value'][line_no])

        # if the values we have are all the same, take the first one and move on

        if min(all_acts) == max(all_acts):
            working_copy = pd.concat([working_copy,sub_lines[0]])
            

        # otherwise look at all the values

        else:
            if min(all_acts) > 1 and max(all_acts) < 11:  # if all the values are "exact values"
                working_copy = pd.concat([working_copy, sub_lines[all_acts.index(min(all_acts))]])  # take the lowest one
            else:
                order = sorted(range(len(all_acts)), key=lambda k: all_acts[k]) # otherwise sort the values
                for item in order:
                    if all_acts[item] != 1:
                        working_copy = pd.concat([working_copy, sub_lines[item]])  # and take the lowest one that isn't 1
                        continue
                        
    in_lines = in_lines.drop(list(set(seen)))
    in_lines = pd.concat([in_lines,working_copy])
    
    if Verbose:
        print('Number of compounds after deduplication pass: ', len(in_lines))

    return in_lines.reset_index(drop=True)




def calc_pscale(value, units):
    if units == 'fM':
        conversion_factor = 10 ** -15
    elif units == 'pM':
        conversion_factor = 10 ** -12
    elif units == 'nM':
        conversion_factor = 10 ** -9
    elif units == 'uM':
        conversion_factor = 10 ** -6
    elif units == 'mM':
        conversion_factor = 10 ** -3
    elif units == 'M':
        conversion_factor = 1
    else:
        print('Unknown units to convert: ', units)
        return -1
    pvalue = (int((math.log10(value * conversion_factor) * -1) * 100)) / 100.0

    return pvalue





def write_smi_act_reg(final_data, base_name, output_dir ='./', add_extra=False):
    nans = final_data[final_data['pchembl_value'].isnull()]
    if len(nans) > 0:
        final_data = final_data.drop(list(nans.index))
        ki = nans[nans['standard_type'] == 'Ki'].reset_index()
        pki = nans[nans['standard_type'].apply(lambda x: any(x in s for s in ['pKi','Log Ki']))].reset_index()
        for i in range(len(ki)):
            ki['pchembl_value'] = calc_pscale(ki['standard_value'],
                                              ki['standard_units'])
        for i in range(len(pki)):
            pki['pchembl_value'] = pki['standard_value']
        final_data = pd.concat([final_data,ki,pki]).reset_index(drop=True)
    
    sort = []
    for i in range(len(final_data['chembl_id'])):
        sort.append(int(final_data['chembl_id'][i][6:]))
    final_data['sort'] = sort
    final_data = final_data.sort_values(by='sort')
    
    smiles_file = base_name + '.smi'
    activity_file = base_name + '.act'
        
    
    df_struct = final_data[['canonical_smiles','chembl_id']]
    df_act = final_data[['chembl_id','pchembl_value']]
    
    #confirm whether or not we want to have the headers written to the files
    df_struct.to_csv(output_dir+'/'+smiles_file, sep = '\t', index = False, header = False)
    df_act.to_csv(output_dir+'/'+activity_file, sep = '\t', index = False, header = False)
    
    return


def write_smi_act_class(buffer, base_name, output_dir, inact_val = 5.0, act_val = 6.0, Verbose = False):
    active_mols = buffer[buffer.pchembl_value >= act_val]
    active_mols['pchembl_value_class'] = np.ones(len(active_mols), dtype=int)
    inactive_mols = buffer[buffer.pchembl_value <= inact_val]
    inactive_mols['pchembl_value_class'] = np.zeros(len(inactive_mols), dtype=int)
    last_list = pd.concat([active_mols,inactive_mols])

    if Verbose:
        print("Total number of molecules in training set: ", len(last_list))
        print("Number of active molecules: ", len(active_mols))
        print("Number of inactive molecules: ", len(inactive_mols))

    df_struct = last_list[['canonical_smiles','chembl_id']]
    df_act = last_list[['chembl_id','pchembl_value_class']]
    df_csv = last_list[['chembl_id','pchembl_value', 'pchembl_value_class','canonical_smiles']]

    smiles_file = base_name + '_class.smi'
    activity_file = base_name + '_class.act'
    output_csv_file = base_name + '_class.csv'
    
    df_struct.to_csv(output_dir+'/'+smiles_file, sep = '\t', index = False, header = False)
    df_act.to_csv(output_dir+'/'+activity_file, sep = '\t', index = False, header = False)
    df_csv.to_csv(output_dir+'/'+output_csv_file, sep = '\t', index = False, header = False)

