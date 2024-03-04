#
#
# filter_config.py
#
# Andrew Fant
#
# 29 March 2019
#
#
# This file holds long constants for the ChEMBL filter programs.  It has been updated to reflect the changes in
# the ChEMBL schema made with release number 25
#

DEBUG = False
VERBOSE = False

Null = ''

kickouts_file_header = ['pref_name', 'organism', 'assay_id', 'assay_type', 'relationship_type', 'relationship_desc',
                        'confidence_score', 'curated_by', 'description', 'activity_id', 'relation', 'value', 'units',
                        'type', 'standard_relation', 'standard_value', 'standard_units', 'standard_flag',
                        'standard_type', 'pchembl_value', 'activity_comment', 'data_validity_comment',
                        'potential_duplicate', 'text_value', 'standard_text_value', 'molregno', 'chembl_id',
                        'canonical_smiles', 'pref_name', 'parent_molregno', 'active_molregno', 'doc_id', 'pubmed_id',
                        'doi', 'journal', 'year', 'volume', 'first_page', 'src_short_name']

test_set_2_compounds = ['CHEMBL633', 'CHEMBL1083993', 'CHEMBL216419', 'CHEMBL1008', 'CHEMBL1200382', 'CHEMBL161',
                        'CHEMBL71', 'CHEMBL1713', 'CHEMBL823', 'CHEMBL799', 'CHEMBL1729', 'CHEMBL1200788',
                        'CHEMBL74656', 'CHEMBL4096162', 'CHEMBL1546244', 'CHEMBL42', 'CHEMBL538973', 'CHEMBL1421',
                        'CHEMBL12', 'CHEMBL543191', 'CHEMBL23', 'CHEMBL1697', 'CHEMBL1200805', 'CHEMBL517',
                        'CHEMBL1201020', 'CHEMBL473', 'CHEMBL502', 'CHEMBL1678', 'CHEMBL1108', 'CHEMBL1175',
                        'CHEMBL1200378', 'CHEMBL652', 'CHEMBL1200822', 'CHEMBL398673', 'CHEMBL1107',
                        'CHEMBL1200901', 'CHEMBL54', 'CHEMBL1200986', 'CHEMBL545608', 'CHEMBL533', 'CHEMBL2355456',
                        'CHEMBL141', 'CHEMBL1230', 'CHEMBL126', 'CHEMBL998', 'CHEMBL651', 'CHEMBL1200825',
                        'CHEMBL137', 'CHEMBL1466172', 'CHEMBL45816', 'CHEMBL1534525', 'CHEMBL58', 'CHEMBL1417019',
                        'CHEMBL32', 'CHEMBL1200735', 'CHEMBL193', 'CHEMBL255863', 'CHEMBL1201740', 'CHEMBL475534',
                        'CHEMBL251230', 'CHEMBL268609', 'CHEMBL1621', 'CHEMBL2107360', 'CHEMBL490', 'CHEMBL1708',
                        'CHEMBL1200609', 'CHEMBL1256912', 'CHEMBL1628650', 'CHEMBL448', 'CHEMBL971', 'CHEMBL16',
                        'CHEMBL553532', 'CHEMBL1611', 'CHEMBL1423', 'CHEMBL702', 'CHEMBL1200820', 'CHEMBL640',
                        'CHEMBL605', 'CHEMBL1294', 'CHEMBL1201486', 'CHEMBL2165709', 'CHEMBL1200437',
                        'CHEMBL3707183', 'CHEMBL254316', 'CHEMBL518520', 'CHEMBL1643', 'CHEMBL1235764',
                        'CHEMBL398463', 'CHEMBL85', 'CHEMBL364512', 'CHEMBL114', 'CHEMBL282042', 'CHEMBL12713',
                        'CHEMBL1422', 'CHEMBL1201174', 'CHEMBL1734', 'CHEMBL1200803', 'CHEMBL471', 'CHEMBL1700',
                        'CHEMBL850', 'CHEMBL535', 'CHEMBL1567', 'CHEMBL374731', 'CHEMBL17157', 'CHEMBL363295',
                        'CHEMBL539770', 'CHEMBL479', 'CHEMBL1200916', 'CHEMBL6966', 'CHEMBL1280', 'CHEMBL638',
                        'CHEMBL2442704']

watch_compounds = ['CHEMBL3', 'CHEMBL4', 'CHEMBL6', 'CHEMBL8', 'CHEMBL22', 'CHEMBL31', 'CHEMBL40', 'CHEMBL79',
                   'CHEMBL108', 'CHEMBL561', 'CHEMBL1087', 'CHEMBL1513', 'CHEMBL1688', 'CHEMBL52348', 'CHEMBL52794',
                   'CHEMBL53321', 'CHEMBL53866', 'CHEMBL284348', 'CHEMBL299390', 'CHEMBL300907', 'CHEMBL412663',
                   'CHEMBL1232472', 'CHEMBL1439973', 'CHEMBL1671895', 'CHEMBL1671896', 'CHEMBL1671899',
                   'CHEMBL1673435', 'CHEMBL2029729', 'CHEMBL2158050', 'CHEMBL3317856']
