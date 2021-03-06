from rdkit.Chem import Descriptors

# Set up descriptor calculation functions and names
dlist = [ Descriptors.BalabanJ,
          Descriptors.BertzCT,
          Descriptors.Chi0,
          Descriptors.Chi0n,
          Descriptors.Chi0v,
          Descriptors.Chi1,
          Descriptors.Chi1n,
          Descriptors.Chi1v,
          Descriptors.Chi2n,
          Descriptors.Chi2v,
          Descriptors.Chi3n,
          Descriptors.Chi3v,
          Descriptors.Chi4n,
          Descriptors.Chi4v,
          Descriptors.EState_VSA1,
          Descriptors.EState_VSA2,
          Descriptors.EState_VSA3,
          Descriptors.EState_VSA4,
          Descriptors.EState_VSA5,
          Descriptors.EState_VSA6,
          Descriptors.EState_VSA7,
          Descriptors.EState_VSA8,
          Descriptors.EState_VSA9,
          Descriptors.EState_VSA10,
          Descriptors.EState_VSA11,
          Descriptors.ExactMolWt,
          Descriptors.FpDensityMorgan1,
          Descriptors.FpDensityMorgan2,
          Descriptors.FpDensityMorgan3,
          Descriptors.FractionCSP3,
          Descriptors.HallKierAlpha,
          Descriptors.HeavyAtomCount,
          Descriptors.HeavyAtomMolWt,
          Descriptors.Ipc,
          Descriptors.Kappa1,
          Descriptors.Kappa2,
          Descriptors.Kappa3,
          Descriptors.LabuteASA,
          Descriptors.MaxAbsEStateIndex,
          Descriptors.MaxAbsPartialCharge,
          Descriptors.MaxEStateIndex,
          Descriptors.MaxPartialCharge,
          Descriptors.MinAbsEStateIndex,
          Descriptors.MinAbsPartialCharge,
          Descriptors.MinEStateIndex,
          Descriptors.MinPartialCharge,
          Descriptors.MolLogP,
          Descriptors.MolMR,
          Descriptors.MolWt,
          Descriptors.NHOHCount,
          Descriptors.NOCount,
          Descriptors.NumAliphaticCarbocycles,
          Descriptors.NumAliphaticHeterocycles,
          Descriptors.NumAliphaticRings,
          Descriptors.NumAromaticCarbocycles,
          Descriptors.NumAromaticHeterocycles,
          Descriptors.NumAromaticRings,
          Descriptors.NumHAcceptors,
          Descriptors.NumHDonors,
          Descriptors.NumHeteroatoms,
          Descriptors.NumRadicalElectrons,
          Descriptors.NumRotatableBonds,
          Descriptors.NumSaturatedCarbocycles,
          Descriptors.NumSaturatedHeterocycles,
          Descriptors.NumSaturatedRings,
          Descriptors.NumValenceElectrons,
          Descriptors.PEOE_VSA1,
          Descriptors.PEOE_VSA2,
          Descriptors.PEOE_VSA3,
          Descriptors.PEOE_VSA4,
          Descriptors.PEOE_VSA5,
          Descriptors.PEOE_VSA6,
          Descriptors.PEOE_VSA7,
          Descriptors.PEOE_VSA8,
          Descriptors.PEOE_VSA9,
          Descriptors.PEOE_VSA10,
          Descriptors.PEOE_VSA11,
          Descriptors.PEOE_VSA12,
          Descriptors.PEOE_VSA13,
          Descriptors.PEOE_VSA14,
          Descriptors.RingCount,
          Descriptors.SMR_VSA1,
          Descriptors.SMR_VSA2,
          Descriptors.SMR_VSA3,
          Descriptors.SMR_VSA4,
          Descriptors.SMR_VSA5,
          Descriptors.SMR_VSA6,
          Descriptors.SMR_VSA7,
          Descriptors.SMR_VSA8,
          Descriptors.SMR_VSA9,
          Descriptors.SMR_VSA10,
          Descriptors.SlogP_VSA1,
          Descriptors.SlogP_VSA2,
          Descriptors.SlogP_VSA3,
          Descriptors.SlogP_VSA4,
          Descriptors.SlogP_VSA5,
          Descriptors.SlogP_VSA6,
          Descriptors.SlogP_VSA7,
          Descriptors.SlogP_VSA8,
          Descriptors.SlogP_VSA9,
          Descriptors.SlogP_VSA10,
          Descriptors.SlogP_VSA11,
          Descriptors.SlogP_VSA12,
          Descriptors.TPSA,
          Descriptors.VSA_EState1,
          Descriptors.VSA_EState2,
          Descriptors.VSA_EState3,
          Descriptors.VSA_EState4,
          Descriptors.VSA_EState5,
          Descriptors.VSA_EState6,
          Descriptors.VSA_EState7,
          Descriptors.VSA_EState8,
          Descriptors.VSA_EState9,
          Descriptors.VSA_EState10,
          Descriptors.fr_Al_COO,
          Descriptors.fr_Al_OH,
          Descriptors.fr_Al_OH_noTert,
          Descriptors.fr_ArN,
          Descriptors.fr_Ar_COO,
          Descriptors.fr_Ar_N,
          Descriptors.fr_Ar_NH,
          Descriptors.fr_Ar_OH,
          Descriptors.fr_COO,
          Descriptors.fr_COO2,
          Descriptors.fr_C_O,
          Descriptors.fr_C_O_noCOO,
          Descriptors.fr_C_S,
          Descriptors.fr_HOCCN,
          Descriptors.fr_Imine,
          Descriptors.fr_NH0,
          Descriptors.fr_NH1,
          Descriptors.fr_NH2,
          Descriptors.fr_N_O,
          Descriptors.fr_Ndealkylation1,
          Descriptors.fr_Ndealkylation2,
          Descriptors.fr_Nhpyrrole,
          Descriptors.fr_SH,
          Descriptors.fr_aldehyde,
          Descriptors.fr_alkyl_carbamate,
          Descriptors.fr_alkyl_halide,
          Descriptors.fr_allylic_oxid,
          Descriptors.fr_amide,
          Descriptors.fr_amidine,
          Descriptors.fr_aniline,
          Descriptors.fr_aryl_methyl,
          Descriptors.fr_azide,
          Descriptors.fr_azo,
          Descriptors.fr_barbitur,
          Descriptors.fr_benzene,
          Descriptors.fr_benzodiazepine,
          Descriptors.fr_bicyclic,
          Descriptors.fr_diazo,
          Descriptors.fr_dihydropyridine,
          Descriptors.fr_epoxide,
          Descriptors.fr_ester,
          Descriptors.fr_ether,
          Descriptors.fr_furan,
          Descriptors.fr_guanido,
          Descriptors.fr_halogen,
          Descriptors.fr_hdrzine,
          Descriptors.fr_hdrzone,
          Descriptors.fr_imidazole,
          Descriptors.fr_imide,
          Descriptors.fr_isocyan,
          Descriptors.fr_isothiocyan,
          Descriptors.fr_ketone,
          Descriptors.fr_ketone_Topliss,
          Descriptors.fr_lactam,
          Descriptors.fr_lactone,
          Descriptors.fr_methoxy,
          Descriptors.fr_morpholine,
          Descriptors.fr_nitrile,
          Descriptors.fr_nitro,
          Descriptors.fr_nitro_arom,
          Descriptors.fr_nitro_arom_nonortho,
          Descriptors.fr_nitroso,
          Descriptors.fr_oxazole,
          Descriptors.fr_oxime,
          Descriptors.fr_para_hydroxylation,
          Descriptors.fr_phenol,
          Descriptors.fr_phenol_noOrthoHbond,
          Descriptors.fr_phos_acid,
          Descriptors.fr_phos_ester,
          Descriptors.fr_piperdine,
          Descriptors.fr_piperzine,
          Descriptors.fr_priamide,
          Descriptors.fr_prisulfonamd,
          Descriptors.fr_pyridine,
          Descriptors.fr_quatN,
          Descriptors.fr_sulfide,
          Descriptors.fr_sulfonamd,
          Descriptors.fr_sulfone,
          Descriptors.fr_term_acetylene,
          Descriptors.fr_tetrazole,
          Descriptors.fr_thiazole,
          Descriptors.fr_thiocyan,
          Descriptors.fr_thiophene,
          Descriptors.fr_unbrch_alkane,
          Descriptors.fr_urea,
          Descriptors.qed
]

dnames = ['BalabanJ',
          'BertzCT',
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
          'EState_VSA1',
          'EState_VSA2',
          'EState_VSA3',
          'EState_VSA4',
          'EState_VSA5',
          'EState_VSA6',
          'EState_VSA7',
          'EState_VSA8',
          'EState_VSA9',
          'EState_VSA10',
          'EState_VSA11',
          'ExactMolWt',
          'FpDensityMorgan1',
          'FpDensityMorgan2',
          'FpDensityMorgan3',
          'FractionCSP3',
          'HallKierAlpha',
          'HeavyAtomCount',
          'HeavyAtomMolWt',
          'Ipc',
          'Kappa1',
          'Kappa2',
          'Kappa3',
          'LabuteASA',
          'MaxAbsEStateIndex',
          'MaxAbsPartialCharge',
          'MaxEStateIndex',
          'MaxPartialCharge',
          'MinAbsEStateIndex',
          'MinAbsPartialCharge',
          'MinEStateIndex',
          'MinPartialCharge',
          'MolLogP',
          'MolMR',
          'MolWt',
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
          'NumRadicalElectrons',
          'NumRotatableBonds',
          'NumSaturatedCarbocycles',
          'NumSaturatedHeterocycles',
          'NumSaturatedRings',
          'NumValenceElectrons',
          'PEOE_VSA1',
          'PEOE_VSA2',
          'PEOE_VSA3',
          'PEOE_VSA4',
          'PEOE_VSA5',
          'PEOE_VSA6',
          'PEOE_VSA7',
          'PEOE_VSA8',
          'PEOE_VSA9',
          'PEOE_VSA10',
          'PEOE_VSA11',
          'PEOE_VSA12',
          'PEOE_VSA13',
          'PEOE_VSA14',
          'RingCount',
          'SMR_VSA1',
          'SMR_VSA2',
          'SMR_VSA3',
          'SMR_VSA4',
          'SMR_VSA5',
          'SMR_VSA6',
          'SMR_VSA7',
          'SMR_VSA8',
          'SMR_VSA9',
          'SMR_VSA10',
          'SlogP_VSA1',
          'SlogP_VSA2',
          'SlogP_VSA3',
          'SlogP_VSA4',
          'SlogP_VSA5',
          'SlogP_VSA6',
          'SlogP_VSA7',
          'SlogP_VSA8',
          'SlogP_VSA9',
          'SlogP_VSA10',
          'SlogP_VSA11',
          'SlogP_VSA12',
          'TPSA',
          'VSA_EState1',
          'VSA_EState2',
          'VSA_EState3',
          'VSA_EState4',
          'VSA_EState5',
          'VSA_EState6',
          'VSA_EState7',
          'VSA_EState8',
          'VSA_EState9',
          'VSA_EState10',
          'fr_Al_COO',
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
          'fr_urea',
          'qed'
]