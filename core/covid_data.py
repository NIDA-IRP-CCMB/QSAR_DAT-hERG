# lsclass=list(df['class'].unique())
lsclass = ['TP', 'FP', 'FN']
color = ['green', 'cyan', 'orange']
# https://matplotlib.org/api/markers_api.html#module-matplotlib.markers
marker = ['o', '+', 'x']
# lsmethodsars=list(df['Method_SARS'].unique())
lsmethodsars = ['PP', 'CPE']
lsfunction = ['antagonist', 'agonist']
simstring = 'sim_to_'

E2_string = '17beta-Estradiol'
tamoxifen_string = 'Tamoxifen'
raloxifene_string = 'Raloxifene'
topotecan_string = 'Topotecan hydrochloride'
# the major androgen
testosterone_string = 'testosterone'

testosterone_smi = 'O=C4\C=C2/[C@]([C@H]1CC[C@@]3([C@@H](O)CC[C@H]3[C@@H]1CC2)C)(C)CC4'

fp_morgan = 'morgan'
fp_rdkit = 'rdkit'
sim_tanimoto = 'tanimoto'
sim_dice = 'dice'


# Descriptor groups are defined below, and config for the groups are defined after

# Group 1: fragments
desc_grp1 = ['fr_Al_COO',
             'fr_Al_OH',
             'fr_Al_OH_noTert',
             'fr_ArN',
             'fr_Ar_COO',
             'fr_Ar_N',
             'fr_Ar_NH',
             'fr_Ar_OH',
             'fr_COO',
             'fr_COO2',
             'fr_C_O',
             'fr_C_O_noCOO',
             'fr_C_S',
             'fr_HOCCN',
             'fr_Imine',
             'fr_NH0',
             'fr_NH1',
             'fr_NH2',
             'fr_N_O',
             'fr_Ndealkylation1',
             'fr_Ndealkylation2',
             'fr_Nhpyrrole',
             'fr_SH',
             'fr_aldehyde',
             'fr_alkyl_carbamate',
             'fr_alkyl_halide',
             'fr_allylic_oxid',
             'fr_amide',
             'fr_amidine',
             'fr_aniline',
             'fr_aryl_methyl',
             'fr_azide',
             'fr_azo',
             'fr_barbitur',
             'fr_benzene',
             'fr_benzodiazepine',
             'fr_bicyclic',
             'fr_diazo',
             'fr_dihydropyridine',
             'fr_epoxide',
             'fr_ester',
             'fr_ether',
             'fr_furan',
             'fr_guanido',
             'fr_halogen',
             'fr_hdrzine',
             'fr_hdrzone',
             'fr_imidazole',
             'fr_imide',
             'fr_isocyan',
             'fr_isothiocyan',
             'fr_ketone',
             'fr_ketone_Topliss',
             'fr_lactam',
             'fr_lactone',
             'fr_methoxy',
             'fr_morpholine',
             'fr_nitrile',
             'fr_nitro',
             'fr_nitro_arom',
             'fr_nitro_arom_nonortho',
             'fr_nitroso',
             'fr_oxazole',
             'fr_oxime',
             'fr_para_hydroxylation',
             'fr_phenol',
             'fr_phenol_noOrthoHbond',
             'fr_phos_acid',
             'fr_phos_ester',
             'fr_piperdine',
             'fr_piperzine',
             'fr_priamide',
             'fr_prisulfonamd',
             'fr_pyridine',
             'fr_quatN',
             'fr_sulfide',
             'fr_sulfonamd',
             'fr_sulfone',
             'fr_term_acetylene',
             'fr_tetrazole',
             'fr_thiazole',
             'fr_thiocyan',
             'fr_thiophene',
             'fr_unbrch_alkane',
             'fr_urea']

# Group 2: Counts
desc_grp2 = ['HeavyAtomCount',
             'NHOHCount',
             'NOCount',
             'NumAliphaticCarbocycles',
             'NumAliphaticHeterocycles',
             'NumAliphaticRings',
             'NumAromaticCarbocycles',
             'NumAromaticHeterocycles',
             'NumAromaticRings',
             'NumHAcceptors',
             'NumHDonors',
             'NumHeteroatoms',
             'NumRotatableBonds',
             'NumSaturatedCarbocycles',
             'NumSaturatedHeterocycles',
             'NumSaturatedRings',
             'RingCount']


# Group 3: Estate and VSA
desc_grp3 = ['EState_VSA1',
             'EState_VSA10',
             'EState_VSA11',
             'EState_VSA2',
             'EState_VSA3',
             'EState_VSA4',
             'EState_VSA5',
             'EState_VSA6',
             'EState_VSA7',
             'EState_VSA8',
             'EState_VSA9',
             'VSA_EState1',
             'VSA_EState10',
             'VSA_EState2',
             'VSA_EState3',
             'VSA_EState4',
             'VSA_EState5',
             'VSA_EState6',
             'VSA_EState7',
             'VSA_EState8',
             'VSA_EState9']

# Group 4: Graph descriptors
desc_grp4 = ['BalabanJ',
             'Chi0',
             'Chi0n',
             'Chi0v',
             'Chi1',
             'Chi1n',
             'Chi1v',
             'Chi2n',
             'Chi2v',
             'Chi3n',
             'Chi3v',
             'Chi4n',
             'Chi4v',
             'Kappa1',
             'Kappa2',
             'Kappa3']

# Group 5: Surface area descriptors
desc_grp5 = ['LabuteASA',
             'PEOE_VSA1',
             'PEOE_VSA10',
             'PEOE_VSA11',
             'PEOE_VSA12',
             'PEOE_VSA13',
             'PEOE_VSA14',
             'PEOE_VSA2',
             'PEOE_VSA3',
             'PEOE_VSA4',
             'PEOE_VSA5',
             'PEOE_VSA6',
             'PEOE_VSA7',
             'PEOE_VSA8',
             'PEOE_VSA9',
             'SMR_VSA1',
             'SMR_VSA10',
             'SMR_VSA2',
             'SMR_VSA3',
             'SMR_VSA4',
             'SMR_VSA5',
             'SMR_VSA6',
             'SMR_VSA7',
             'SMR_VSA8',
             'SMR_VSA9',
             'SlogP_VSA1',
             'SlogP_VSA10',
             'SlogP_VSA11',
             'SlogP_VSA12',
             'SlogP_VSA2',
             'SlogP_VSA3',
             'SlogP_VSA4',
             'SlogP_VSA5',
             'SlogP_VSA6',
             'SlogP_VSA7',
             'SlogP_VSA8',
             'SlogP_VSA9',
             'TPSA']

# Group 6: Decimals and ratios
desc_grp6 = ['FractionCSP3',
             'qed']

# Group 7: Weirdos
desc_grp7 = ['BertzCT',
             'HallKierAlpha',
             'Ipc']

