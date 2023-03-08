
ribosome_stoich = {'30_S_assembly':{'stoich':{'generic_16s_rRNAs':1}},
                    '50_S_assembly':{'stoich':{'generic_23s_rRNAs':1,'generic_5s_rRNAs':1}},
                    'assemble_ribosome_subunits': {'stoich': {'gtp_c': 1}}}

ribosome_subreactions = {
    'gtp_bound_30S_assembly_factor_phase1':
                         {'enzyme': '',
                          'stoich': {"gtp_c": -2,
                                     "h2o_c": -2,
                                     "h_c": 2,
                                     "pi_c": 2,
                                     "gdp_c": 2},
                          'num_mods': 1},
                         'RbfA_mono_assembly_factor_phase1':
                         {'enzyme': '',
                          'stoich': {},
                          'num_mods': 1},
                         'RimM_mono_assembly_factor_phase1':
                         {'enzyme': '',
                          'stoich': {},
                          'num_mods': 1}
}

generics = {'generic_16Sm4Cm1402': {'enzymes':[]},
                        'generic_LYSINEaaRS': {'enzymes':[]},
                        'generic_Dus': {'enzymes':[]},
                        'generic_RF': {'enzymes':[]},
                        'generic_Tuf': {'enzymes':[]},
                        'generic_RNase': {'enzymes':[]},
                        'generic_16s_rRNAs': {'enzymes':[]},
                        'generic_23s_rRNAs': {'enzymes':[]},
                        'generic_5s_rRNAs': {'enzymes':[]},
                        'generic_5s_rRNAs': {'enzymes':[]},
                        'generic_fes_transfers_complex': {'enzymes':[]}
           }

rrna_modifications = {
                      # ---------16S Modifications---------------
                      'm2G_at_1207': {'machine': '',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm2G_at_1516': {'machine': '',
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm2G_at_966': {'machine': '',  
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm3U_at_1498': {'machine': '',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm4Cm_at_1402': {'machine': '',
                                       'metabolites': {'amet_c': -2,
                                                       'ahcys_c': 2,
                                                       'h_c': 2}},
                      'm5C_at_1407': {'machine': '',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5C_at_967': {'machine': '',  
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm62A_at_1518': {'machine': '',  
                                       'metabolites': {'amet_c': -2,
                                                       'ahcys_c': 2,
                                                       'h_c': 2}},
                      'm62A_at_1519': {'machine': '',  
                                       'metabolites': {'amet_c': -2,
                                                       'ahcys_c': 2,
                                                       'h_c': 2}},
                      'm7G_at_527': {'machine': '',  
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'Y_at_516': {'machine': '',  
                                   'metabolites': {}},

                      # ---------23S Modifications---------------
                      'Cm_at_2498': {'machine': '',
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'D_at_2449': {'machine': '',
                                    'metabolites': {'h_c': -1,
                                                    'nadh_c': -1,
                                                    'nad_c': 1}},
                      'Gm_at_2251': {'machine': '',  
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm1G_at_745': {'machine': '',
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm2A_at_2503': {'machine': '',
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm2G_at_1835': {'machine': '',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm2G_at_2445': {'machine': '',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5C_at_1962': {'machine': '',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5U_at_1939': {'machine': '',
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5U_at_747': {'machine': '',
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm6A_at_1618': {'machine': '',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm6A_at_2030': {'machine': '',
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm7G_at_2069': {'machine': '',
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'Um_at_2552': {'machine': '',  
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'Y_at_1911': {'machine': '',
                                    'metabolites': {}},
                      'Y_at_1915': {'machine': '',
                                    'metabolites': {}},
                      'm3Y_at_1915': {'machine': '',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'Y_at_1917': {'machine': '',
                                    'metabolites': {}},
                      'Y_at_2457': {'machine': '',  
                                    'metabolites': {}},
                      'Y_at_2504': {'machine': '',  
                                    'metabolites': {}},
                      'Y_at_2580': {'machine': '',  
                                    'metabolites': {}},
                      'Y_at_2604': {'machine': '',  
                                    'metabolites': {}},
                      'Y_at_2605': {'machine': '',  
                                    'metabolites': {}},
                      'Y_at_746': {'machine': '',  
                                   'metabolites': {}},
                      'Y_at_955': {'machine': '',  
                                   'metabolites': {}}}

rrna_modification_info = {'Y': {'elements': {}, 'charge': 0},
                     'm3Y': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'Um': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'm7G': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm6A': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm5U': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm5C': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm2G': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm2A': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm1G': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'Gm': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'D': {'elements': {'H': 2}, 'charge': 0},
                     'Cm': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm62A': {'elements': {'C': 2, 'H': 4}, 'charge': 0},
                     'm4Cm': {'elements': {'C': 2, 'H': 4}, 'charge': 0},
                     'm3U': {'elements': {'C': 1, 'H': 2}, 'charge': 0}
                     }

amino_acid_trna_synthetase = {
  "cys__L_c": '',
  "leu__L_c": '',
  "lys__L_c": '',
  "asp__L_c": '',
  "phe__L_c": '',
  "his__L_c": '',
  "asn__L_c": '',
  "pro__L_c": '',
  "ala__L_c": '',
  "ile__L_c": '',
  "ser__L_c": '',
  "arg__L_c": '',
  "met__L_c": '',
  "tyr__L_c": '',
  "glu__L_c": '',
  "thr__L_c": '',
  "val__L_c": '',
  "gly_c": '',
  "trp__L_c": '',
  "gln__L_c": '',
}

amino_acid_regex = {'cys__L_c': 'cystein',
 'leu__L_c': '(?<!iso)leuc',
 'lys__L_c': 'lys',
 'asp__L_c': 'aspart',
 'phe__L_c': 'phenylalan',
 'his__L_c': 'histid',
 'asn__L_c': 'asparag',
 'pro__L_c': 'prol',
 'ala__L_c': 'alan',
 'ile__L_c': 'isoleuc',
 'ser__L_c': 'ser',
 'arg__L_c': 'argin',
 'met__L_c': 'methion',
 'tyr__L_c': 'tyros',
 'glu__L_c': 'glutam(?!inyl|ine)',
 'thr__L_c': 'threon',
 'val__L_c': 'val',
 'gly_c': 'glyc',
 'trp__L_c': 'tryptophan',
 'gln__L_c': 'glutamin'}

translation_stop_dict = {'UAG': {'enzyme':''},
                         'UGA': {'enzyme':''},
                         'UAA': {'enzyme':'generic_RF'}}

initiation_subreactions = {
    'Translation_initiation_factor_InfA':
        {'enzymes': [],
         'stoich': {},
        'element_contribution': {}},

    'Translation_initiation_factor_InfC':
        {'enzymes': [],
         'stoich': {},
        'element_contribution': {}},

    'Translation_gtp_initiation_factor_InfB':
        {'enzymes': [],
         'stoich': {'gtp_c': -1,
                    'h2o_c': -1,
                    'h_c': 1,
                    'pi_c': 1,
                    'gdp_c': 1},
        'element_contribution': {}},

    'fmet_addition_at_START':
        {'enzymes': [],
         # iOL had h_c:1 for fmet addition but this is not mass balanced
         'stoich': {'10fthf_c': -1, 'thf_c': 1,
                    # 'h_c': 1,
                    'generic_tRNA_START_met__L_c': -1},
         'element_contribution': {'C': 1, 'O': 1}}
   }

elongation_subreactions = {'FusA_mono_elongation': {'enzymes': [],
                                                    'stoich': {'gtp_c': -1,
                                                               'h2o_c': -1,
                                                               'h_c': 1,
                                                               'pi_c': 1,
                                                               'gdp_c': 1}},

                           'Tuf_gtp_regeneration': {'enzymes': [],
                                                    'stoich': {}}}
termination_subreactions = {'PrfA_mono_mediated_termination':
                            {'enzymes': [],
                             'stoich': {}},

                            'PrfB_mono_mediated_termination':
                            {'enzymes': [],
                             'stoich': {}},

                            'generic_RF_mediated_termination':
                            {'enzymes': [],
                             'stoich': {}},

                            'N_terminal_methionine_cleavage':
                            {'enzymes': [],
                             'stoich': {'h2o_c': -1,
                                        'met__L_c': 1, 'h_c': 1},
                             'element_contribution': {'H': -10, 'O': -1,
                                                      'C': -5, 'N': -1,
                                                      'S': -1}},

                            'peptide_deformylase_processing':
                            {'enzymes': [],
                             'stoich': {'h2o_c': -1,
                                        'for_c': 1},
                             'element_contribution':
                                 {'H': 1, 'O': -1, 'C': -1}},

                            # This is a GTPS
                            'peptide_chain_release':
                            {'enzymes': [],
                             'stoich': {'gtp_c': -1,
                                        'h2o_c': -1,
                                        'h_c': 1,
                                        'pi_c': 1,
                                        'gdp_c': 1}},

                            'ribosome_recycler':
                            {'enzymes': [],
                             'stoich': {}},

                            'GroEL_dependent_folding':
                            {'enzymes': [],
                             'stoich': {'atp_c': -7,
                                        'h2o_c': -7,
                                        'h_c': 7,
                                        'adp_c': 7,
                                        'pi_c': 7}},

                            # DnaK is correct
                            'DnaK_dependent_folding':
                            {'enzymes': [],
                             'stoich': {'atp_c': -1,
                                        'h2o_c': -1,
                                        'h_c': 1,
                                        'adp_c': 1,
                                        'pi_c': 1}}
                            }
special_trna_subreactions = {
    'sec_addition_at_UGA': {
        'enzymes': [],  # Selenocysteine loaders
        'stoich': {'h_c': 1, 'selnp_c': -1,
                   'pi_c': 1,
                   'generic_tRNA_UGA_ser__L_c': -1}, # Serine is the precursor to Sec
        'element_contribution': {'O': -1, 'Se': 1}}}

excision_machinery = {
    'rRNA_containing': {'enzymes':[]},
    'monocistronic': {'enzymes':[]},
    'polycistronic_wout_rRNA': {'enzymes':[]}}

folding_dict = {
    "GroEL_dependent_folding": {'enzymes':[]},
    "DnaK_dependent_folding": {'enzymes':[]}
}

special_modifications = {
    'fes_transfer':{'enzymes':[],'stoich':{}},
    'fes_chaperones':{'enzymes':[],'stoich':{}},
    'generic_2fe2s_transfer_complex':{'enzymes':[],'stoich':{}},
    'generic_4fe4s_transfer_complex':{'enzymes':[],'stoich':{}},
    'lipoate_modifications':{'enzymes':[],'stoich':{}},
    'bmocogdp_chaperones':{'enzymes':[],'stoich':{}}
}

trna_modification = {'D_at_20A': {'enzymes': [],
                                  'stoich': {'nadph_c': -1,
                                                  'h_c': -1,
                                                  'nadp_c': 1}},

                     'D_at_20': {'enzymes': [],  
                                 'stoich': {'nadph_c': -1,
                                                 'h_c': -1,
                                                 'nadp_c': 1}},

                     't6A_at_37': {'enzymes': [],  
                                   'stoich': {'hco3_c': -1,
                                                   'thr__L_c': -1,
                                                   'atp_c': -1,
                                                   'amp_c': 1,
                                                   'h_c': 1,
                                                   'h2o_c': 1,
                                                   'ppi_c': 1}},

                     'm7G_at_46': {'enzymes': [],  
                                   'stoich': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     'acp3U_at_47': {'enzymes': [],
                                     'stoich': {'amet_c': -1,
                                                     '5mta_c': 1,
                                                     'h_c': 1}},

                     'm5U_at_54': {'enzymes': [],  
                                   'stoich': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     'Y_at_55': {'enzymes': [],  
                                 'stoich': {}},

                     'Y_at_65': {'enzymes': [],  
                                 'stoich': {}},

                     'D_at_17': {'enzymes': [],  
                                 'stoich': {'nadph_c': -1,
                                                 'h_c': -1,
                                                 'nadp_c': 1}},

                     'cmo5U_at_34': {'enzymes': [],
                                     'stoich': {'amet_c': -1,
                                                     'chor_c': -2,
                                                     'ahcys_c': 1,
                                                     'h_c': 1,
                                                     'C10H8O5_c': 1,
                                                     'C9H9O4_c': 1}},

                     'D_at_16': {'enzymes': [],  
                                 'stoich': {'nadph_c': -1,
                                                 'h_c': -1,
                                                 'nadp_c': 1}},

                     'Q_at_34': {
                     'enzymes': [],
                     'stoich': {'preq1_c': -1,
                                     'amet_c': -1,
                                     'gua_c': 1,
                                     'ade_c': 1,
                                     'met__L_c': 1,
                                     'h_c': 2}},

                     'm2A_at_37': {'enzymes': [],
                                   'stoich': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     's4U_at_8': {'enzymes': [],
                                  'carriers': {},
                                  'stoich': {'atp_c': -1,
                                                  'amp_c': 1,
                                                  'ppi_c': 1,
                                                  'h_c': 1}},

                     'm6t6A_at_37': {'enzymes': [],
                                     'stoich': {'amet_c': -1,
                                                     'atp_c': -1,
                                                     'hco3_c': -1,
                                                     'thr__L_c': -1,
                                                     'ahcys_c': 1,
                                                     'amp_c': 1,
                                                     'h_c': 2,
                                                     'h2o_c': 1,
                                                     'ppi_c': 1}},
                     's2C_at_32': {'enzymes': [],
                                   'carriers': {},
                                   'stoich': {'atp_c': -1,
                                                   'amp_c': 1,
                                                   'ppi_c': 1,
                                                   'h_c': 1}},
                     'mnm5U_at_34': {'enzymes': [],
                                     'stoich': {'gtp_c': -1,
                                                     'h2o_c': -1,
                                                     '5fthf_c': -1,
                                                     'gly_c': -1,
                                                     'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 3,
                                                     'gdp_c': 1,
                                                     'glx_c': 1,
                                                     'pi_c': 1,
                                                     'thf_c': 1}},

                     'Y_at_40': {'enzymes': [],  
                                 'stoich': {}},

                     'Gm_at_18': {'enzymes': [],  
                                  'stoich': {'amet_c': -1,
                                                  'ahcys_c': 1,
                                                  'h_c': 1}},

                     'Um_at_32': {'enzymes': [],
                                  'stoich': {'amet_c': -1,
                                                  'ahcys_c': 1,
                                                  'h_c': 1}},

                     'Y_at_38': {'enzymes': [],  
                                 'stoich': {}},

                     'ac4C_at_34': {'enzymes': [],  
                                    'stoich': {'accoa_c': -1,
                                                    'coa_c': 1}},

                     'Y_at_39': {'enzymes': [],  
                                 'stoich': {}},
                     
                     'mnm5s2U_at_34': {'enzymes': [],
                                       'carriers': {},
                                       'stoich': {'atp_c': -1,
                                                       'gtp_c': -1,
                                                       'h2o_c': -1,
                                                       '5fthf_c': -1,
                                                       'gly_c': -1,
                                                       'amet_c': -1,
                                                       'gdp_c': 1,
                                                       'pi_c': 1,
                                                       'h_c': 4,
                                                       'thf_c': 1,
                                                       'glx_c': 1,
                                                       'ahcys_c': 1,
                                                       'amp_c': 1,
                                                       'ppi_c': 1}},

                     'm6A_at_37': {'enzymes': [],
                                   'stoich': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     'Cm_at_32': {'enzymes': [],
                                  'stoich': {'amet_c': -1,
                                                  'ahcys_c': 1,
                                                  'h_c': 1}},

                     'ms2i6A_at_37': {'enzymes': [],

                                      'carriers': {},
                                      'stoich': {'dmpp_c': -1,
                                                      'amet_c': -2,
                                                      'ppi_c': 1,
                                                      'ahcys_c': 1,
                                                      'h_c': 2,
                                                      'met__L_c': 1,
                                                      'dad__5_c': 1,
                                                      }},

                     'Y_at_32': {'enzymes': [],  
                                 'stoich': {}},

                     'D_at_21': {'enzymes': [],  
                                 'stoich': {'nadph_c': -1,
                                                 'h_c': -1,
                                                 'nadp_c': 1}},

                     'm1G_at_37': {'enzymes': [],  
                                   'stoich': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     'Y_at_13': {'enzymes': [],  
                                 'stoich': {}},

                     'k2C_at_34': {'enzymes': [],  
                                   'stoich': {'atp_c': -1,
                                                   'lys__L_c': -1,
                                                   'ppi_c': 1,
                                                   'amp_c': 1,
                                                   'h_c': 2}},

                     'I_at_34': {'enzymes': [],
                                 'stoich': {'h2o_c': -1, 'h_c': -1,
                                                 'nh4_c': 1}},

                     'i6A_at_37': {'enzymes': [],  
                                   'stoich': {'dmpp_c': -1,
                                                   'ppi_c': 1}},

                     'D_at_20_in_met_tRNA': {'enzymes': [],
                                             
                                             'stoich': {'nadph_c': -1,
                                                             'h_c': -1,
                                                             'nadp_c': 1}},

                     'D_at_16_in_met_tRNA': {'enzymes': [],
                                             
                                             'stoich': {'nadph_c': -1,
                                                             'h_c': -1,
                                                             'nadp_c': 1}},

                     'D_at_17_in_met_tRNA': {'enzymes': [],
                                             
                                             'stoich': {'nadph_c': -1,
                                                             'h_c': -1,
                                                             'nadp_c': 1}},
                     'D_at_20A_in_met_tRNA': {'enzymes': [],
                                              
                                              'stoich': {'nadph_c': -1,
                                                              'h_c': -1,
                                                              'nadp_c': 1}}
                     }

trna_modification_info = {'D': {'elements': {'H': 2}, 'charge': 0},
                     'i6A': {'elements': {'C': 5, 'H': 8}, 'charge': 0},
                     'I': {'elements': {'N': -1, 'H': -1, 'O': 1},
                           'charge': 0},
                     'k2C': {'elements': {'O': 1, 'N': 2, 'H': 12, 'C': 6},
                             'charge': 0},
                     'Y': {'elements': {}, 'charge': 0},
                     'm1G': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'ms2i6A': {'elements': {'C': 6, 'H': 10, 'S': 1},
                                'charge': 0},
                     'Cm': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'Um': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'm6A': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'mnm5s2U': {'elements': {'C': 2, 'H': 5, 'N': 1, 'O': -1,
                                              'S': 1},
                                 'charge': 0},
                     'ac4C': {'elements': {'H': 2, 'C': 2, 'O': 1},
                              'charge': 0},
                     'Gm': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'mnm5U': {'elements': {'C': 2, 'H': 5, 'N': 1},
                               'charge': 0},
                     's2C': {'elements': {'O': -1, 'S': 1}, 'charge': 0},
                     'm6t6A': {'elements': {'C': 6, 'O': 4, 'N': 1, 'H': 9},
                               'charge': 0},
                     's4U': {'elements': {'O': -1, 'S': 1}, 'charge': 0},
                     'm2A': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'Q': {'elements': {'C': 7, 'O': 2, 'H': 11}, 'charge': 1},
                     'cmo5U': {'elements': {'C': 2, 'O': 3, 'H': 2},
                               'charge': 0},
                     'm5U': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'acp3U': {'elements': {'C': 4, 'H': 7, 'N': 1, 'O': 2},
                               'charge': 0},
                     'm7G': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     't6A': {'elements': {'C': 5, 'N': 1, 'O': 4, 'H': 6},
                             'charge': 0}
                     }

peptide_processing_subreactions = {"peptide_deformylase_processing",
                                   "peptide_chain_release",
                                   "ribosome_recycler"}

translation_start_codons = {"AUG", "GUG", "UUG", "AUU", "CUG"}

transcription_subreactions = {
    'Transcription_normal_rho_independent':
        {'enzymes': [],
         'stoich': {}},
    'Transcription_normal_rho_dependent':
        {'enzymes': [],
         'stoich': {'atp_c': -3,
                    'h2o_c': -3,
                    'adp_c': 3,
                    'pi_c': 3,
                    'h_c': 3}},
    'Transcription_stable_rho_independent':
        {'enzymes': [],
         'stoich': {}},
    'Transcription_stable_rho_dependent':
        {'enzymes': [],
         'stoich': {'atp_c': -3,
                    'h2o_c': -3,
                    'adp_c': 3,
                    'pi_c': 3,
                    'h_c': 3}}
}

abbreviation_to_pathway = {'s': 'sec_translocation',
                           't': ['tat_translocation', 'tat_translocation_alt'],
                           'b': 'bam_translocation',
                           'l': 'lol_translocation',
                           'y': 'yidC_translocation',
                           'a': 'secA_translocation',
                           'p': 'srp_yidC_translocation',
                           'r': 'srp_translocation'}

pathway_to_abbreviation = {}
for a,p in abbreviation_to_pathway.items():
    if isinstance(p,list):
        for i in p:
            pathway_to_abbreviation[i.replace('_translocation','')] = a
    else:
        pathway_to_abbreviation[p.replace('_translocation','')] = a