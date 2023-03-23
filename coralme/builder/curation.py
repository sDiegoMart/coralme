import pandas
import os
import coralme
from coralme.builder import dictionaries
import re

import logging
log = logging.getLogger(__name__)

import tqdm
import cobra

bar_format = '{desc:<75}: {percentage:.1f}%|{bar}| {n_fmt:>5}/{total_fmt:>5} [{elapsed}<{remaining}]'
class MEManualCuration(object):
    """MEManualCuration class for loading manual curation from files

    This class acts as a loader for all manual curation inputs. It is
    used by Organism to retrieve manual curation in a coralME-ready
    format.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.
    """

    def __init__(self,
                 org):
        self.org = org
        self.directory = self.org.directory
        self.curation_notes = self.org.curation_notes
        self.is_reference = self.org.is_reference
        self.configuration = self.org.config

    def load_manual_curation(self):
        logging.warning("Loading protein location")
        self.org.protein_location = self.load_protein_location()
        logging.warning("Loading protein translocation multipliers")
        self.org.translocation_multipliers = self.load_translocation_multipliers()
        logging.warning("Loading lipoprotein precursors")
        self.org.lipoprotein_precursors = self.load_lipoprotein_precursors()
        logging.warning("Loading methionine cleaved proteins")
        self.org.cleaved_methionine = self.load_cleaved_methionine()
        logging.warning("Loading subsystem classification for Keffs")
        self.org.subsystem_classification = self.load_subsystem_classification()
        logging.warning("Loading manually added complexes")
        self.org.manual_complexes = self.load_manual_complexes()
        logging.warning("Loading reaction corrections")
        self.org.reaction_corrections = self.load_reaction_corrections()
        logging.warning("Loading sigma factors")
        self.org.sigmas = self.load_sigmas()
        logging.warning("Loading M to ME metabolites dictionary")
        self.org.me_mets = self.load_me_mets()
        logging.warning("Loading RNA degradosome")
        self.org.rna_degradosome = self.load_rna_degradosome()
        logging.warning("Reading ribosomal proteins")
        self.org.ribosome_stoich = self.load_ribosome_stoich()
        logging.warning("Loading ribosome subreactions")
        self.org.ribosome_subreactions = self.load_ribosome_subreactions()
        logging.warning("Loading generics")
        self.org.generic_dict = self.load_generic_dict()
#         logging.warning("Loading ribosome rrna modifications")
#         self.org.rrna_modifications = self.load_rrna_modifications()
        logging.warning("Loading amino acid tRNA synthetases")
        self.org.amino_acid_trna_synthetase = self.load_amino_acid_trna_synthetase()
        logging.warning("Loading peptide release factors")
        self.org.peptide_release_factors = self.load_peptide_release_factors()
        logging.warning("Loading translation initiation subreactions")
        self.org.initiation_subreactions = self.load_initiation_subreactions()
        logging.warning("Loading translation elongation subreactions")
        self.org.elongation_subreactions = self.load_elongation_subreactions()
        logging.warning("Loading translation termination subreactions")
        self.org.termination_subreactions = self.load_termination_subreactions()
        logging.warning("Loading special trna subreactions")
        self.org.special_trna_subreactions = self.load_special_trna_subreactions()
        logging.warning("Loading RNA excision machinery")
        self.org.excision_machinery = self.load_excision_machinery()
        logging.warning("Loading special modifications")
        self.org.special_modifications = self.load_special_modifications()
        logging.warning("Loading rna modifications and targets")
        self.org.rna_modification = self.load_rna_modification()
        self.org.rna_modification_targets = self.load_rna_modification_targets()
        logging.warning("Loading folding information of proteins")
        self.org.folding_dict = self.load_folding_dict()
        logging.warning("Loading transcription subreactions")
        self.org.transcription_subreactions = self.load_transcription_subreactions()
        logging.warning("Loading protein translocation pathways")
        self.org.translocation_pathways = self.load_translocation_pathways()
        logging.warning("Loading lipid modifications")
        self.org.lipid_modifications = self.load_lipid_modifications()

        # Special input files
        logging.warning("Loading subreaction matrix")
        self.org.subreaction_matrix = self.load_subreaction_matrix()
        logging.warning("Loading reaction matrix")
        self.org.reaction_matrix = self.load_reaction_matrix()
        logging.warning("Loading stable RNAs")
        self.org.stable_RNAs = self.load_stable_RNAs()
        logging.warning("Loading orphan reactions")
        self.org.orphan_and_spont_reactions = self.load_orphan_and_spont_reactions()
        logging.warning("Loading enzyme-reaction-association")
        self.org.enz_rxn_assoc_df = self.load_enz_rxn_assoc_df()


    def _get_manual_curation(self,
                             filename,
                             create_file=None,
                             no_file_return=pandas.DataFrame(),
                             sep = '\t',
                             pathtype = 'relative'):
        if pathtype == 'absolute':
            filepath = filename
        else:
            filepath = self.directory + filename
        if os.path.isfile(filepath):
            return pandas.read_csv(filepath,
                                   index_col=0,
                                   sep=sep,
                                  comment='#',
                                  skip_blank_lines=True).fillna("")

        if create_file is not None:
            create_file.to_csv(filepath,sep=sep)

        self.curation_notes['org._get_manual_curation'].append({
            'msg':'No {} file found'.format(filename),
            'importance':'low',
            'to_do':'Fill in {}'.format(filepath)
        })
        return no_file_return

    def load_reaction_corrections(self):
        create_file = pandas.DataFrame(columns = [
                'reaction_id',
                'name',
                'gene_reaction_rule',
                'reaction',
                'notes',
                ]).set_index('reaction_id')
        return self._get_manual_curation(
             "reaction_corrections.txt",
             create_file = create_file,
             no_file_return = create_file,
             sep = ',').T.to_dict()

    def load_protein_location(self):
        create_file = pandas.DataFrame(columns = [
                'Complex',
                'Complex_compartment',
                'Protein',
                'Protein_compartment',
                'translocase_pathway',
                ]).set_index('Complex')
        return self._get_manual_curation(
             "peptide_compartment_and_pathways.txt",
             create_file = create_file,
             no_file_return = create_file)

    def load_translocation_multipliers(self):
        create_file = None
        return self._get_manual_curation(
             "translocation_multipliers.txt",
            create_file=create_file,
            no_file_return=pandas.DataFrame(),
            sep=',').to_dict()

    def load_lipoprotein_precursors(self):
        create_file = pandas.DataFrame(columns=['gene'])
        return self._get_manual_curation(
            "lipoprotein_precursors.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep=',').to_dict()["gene"]

    def load_cleaved_methionine(self):
        create_file = pandas.DataFrame.from_dict({'cleaved_methionine_genes':{}}).set_index('cleaved_methionine_genes')
        return self._get_manual_curation(
            "cleaved_methionine.txt",
            create_file = create_file,
            no_file_return = create_file).index.to_list()

    def _create_subsystem_classification(self,subsystems):
        d = {}
        for s in subsystems:
            d[s] = {}
            d[s]["central_CE"] = 0
            d[s]["central_AFN"] = 0
            d[s]["intermediate"] = 0
            d[s]["secondary"] = 0
            d[s]["other"] = 1
        return pandas.DataFrame.from_dict(d).T
    def load_subsystem_classification(self):
        if self.is_reference:
            return None

        subsystems = set(r.subsystem for r in self.org.m_model.reactions if r.subsystem)
        create_file = self._create_subsystem_classification(subsystems)

        df = self._get_manual_curation(
            "subsystem_classification.txt",
            create_file = create_file,
            no_file_return = pandas.DataFrame())

        d = {}
        for c in df.columns:
            for s in df[df[c] == 1].index:
                d[s] = c
        return d

    def load_manual_complexes(self):
        create_file = pandas.DataFrame.from_dict(
                {"complex_id": {}, "name": {}, "genes": {}, "mod": {}, "replace": {}}
            ).set_index("complex_id")
        return self._get_manual_curation(
            "protein_corrections.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = ',')

    def load_sigmas(self):
        create_file = pandas.DataFrame(columns = [
                'sigma','complex','genes','name'
            ]).set_index('sigma')
        return self._get_manual_curation(
            "sigma_factors.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = ',')

    def load_me_mets(self):
        create_file = pandas.DataFrame(columns = [
                'id','me_id','name','formula','compartment','type'
            ]).set_index('id')
#         filename = self.configuration.get('df_metadata_metabolites', None)
#         if filename is None or not filename:
#             filename = self.org.directory + "me_metabolites.txt"
        filename = self.org.directory + "me_metabolites.txt"
        return self._get_manual_curation(
            filename,
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t',
            pathtype = 'absolute')

    def load_rna_degradosome(self):
        create_file = pandas.DataFrame(columns = [
                'rna_degradosome'
            ]).set_index('rna_degradosome')
        df = self._get_manual_curation(
            "rna_degradosome.txt",
            create_file = create_file,
            no_file_return = create_file)
        return {
            'rna_degradosome' : {
                'enzymes' : list(df.index)
            }
        }

    def _process_ribosome_stoich(self,
                               df):
        from copy import deepcopy
        ribosome_stoich = deepcopy(dictionaries.ribosome_stoich)
        for s, row in df.iterrows():
            proteins = row["proteins"]
            if proteins:
                proteins = proteins.split(",")
                for p in proteins:
                    if s == "30S":
                        ribosome_stoich["30_S_assembly"]["stoich"][p] = 1
                    elif s == "50S":
                        ribosome_stoich["50_S_assembly"]["stoich"][p] = 1
        return ribosome_stoich

    def _create_ribosome_stoich(self):
        return pandas.DataFrame.from_dict(
            {"proteins": {"30S": "generic_16s_rRNAs",
                          "50S": "generic_5s_rRNAs,generic_23s_rRNAs"}}).rename_axis('subunit')

    def load_ribosome_stoich(self):
        create_file = self._create_ribosome_stoich()
        df = self._get_manual_curation(
            "ribosomal_proteins.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')

        return self._process_ribosome_stoich(df)

    def _create_ribosome_subreactions(self):
        return pandas.DataFrame.from_dict(dictionaries.ribosome_subreactions.copy()).T.rename_axis('subreaction')
    def _modify_ribosome_subreactions_for_save(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "stoich"] = self._dict_to_str(row["stoich"])
        return df
    def _modify_ribosome_subreactions_from_load(self,df):
        d = df.T.to_dict()
        for r, row in d.items():
            d[r]["stoich"] = self._str_to_dict(row["stoich"])
        return d
    def load_ribosome_subreactions(self):
        create_file = self._modify_ribosome_subreactions_for_save(
                self._create_ribosome_subreactions())
        df = self._get_manual_curation(
            "ribosome_subreactions.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_ribosome_subreactions_from_load(df)

    def _create_generic_dict(self):
        return pandas.DataFrame.from_dict(dictionaries.generics.copy()).T.rename_axis('generic_component')
    def _modify_generic_dict_for_save(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_generic_dict_from_load(self,df):
        d = df.T.to_dict()
        for k, v in d.items():
            d[k] = {}
            if v['enzymes']:
                d[k]['enzymes'] = v['enzymes'].split(",")
            else:
                d[k]['enzymes'] = []
        return d
    def load_generic_dict(self):
        create_file = self._modify_generic_dict_for_save(
                self._create_generic_dict())
        df = self._get_manual_curation(
            "generic_dict.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_generic_dict_from_load(df)

#     def _create_rrna_modifications(self):
#         return pandas.DataFrame.from_dict(dictionaries.rrna_modifications.copy()).T.rename_axis('modification')
#     def _modify_rrna_modifications_for_save(self,df):
#         df = df.copy()
#         for r, row in df.iterrows():
#             df.loc[r, "metabolites"] = self._dict_to_str(row["metabolites"])
#         return df
#     def _modify_rrna_modifications_from_load(self,df):
#         d = df.T.to_dict()
#         for k, v in d.items():
#             v["metabolites"] = self._str_to_dict(v["metabolites"])
#         return d
#     def load_rrna_modifications(self):
#         create_file = self._modify_rrna_modifications_for_save(
#                 self._create_rrna_modifications())
#         df = self._get_manual_curation(
#             "rrna_modifications.txt",
#             create_file = create_file,
#             no_file_return = create_file,
#             sep = '\t')
#         return self._modify_rrna_modifications_from_load(df)

    def _create_amino_acid_trna_synthetase(self):
        return pandas.DataFrame.from_dict(
            {"enzyme": dictionaries.amino_acid_trna_synthetase.copy()}).rename_axis('amino_acid')
    def load_amino_acid_trna_synthetase(self):
        create_file = self._create_amino_acid_trna_synthetase()
        return self._get_manual_curation(
                "amino_acid_trna_synthetase.txt",
                create_file = create_file,
                no_file_return = create_file,
                sep = '\t').to_dict()['enzyme']

    def _create_peptide_release_factors(self):
        return pandas.DataFrame.from_dict(dictionaries.translation_stop_dict.copy()).T.rename_axis('release_factor')
    def _modify_peptide_release_factors_for_save(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzyme"] = row["enzyme"]
        return df
    def _modify_peptide_release_factors_from_load(self,df):
        d = df.T.to_dict()
        for k, v in d.items():
            d[k] = {}
            if v['enzyme']:
                d[k]['enzyme'] = v['enzyme']
            else:
                d[k]['enzyme'] = ''
        return d
    def load_peptide_release_factors(self):
        create_file = self._modify_peptide_release_factors_for_save(
                self._create_peptide_release_factors())
        df = self._get_manual_curation(
            "peptide_release_factors.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_peptide_release_factors_from_load(df)

    def _create_initiation_subreactions(self):
        return pandas.DataFrame.from_dict(dictionaries.initiation_subreactions.copy()).T.rename_axis('subreaction')
    def _modify_initiation_subreactions_for_save(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
            df.loc[r, "stoich"] = self._dict_to_str(row["stoich"])
            df.loc[r, "element_contribution"] = self._dict_to_str(
                row["element_contribution"]
            )
        return df
    def _modify_initiation_subreactions_from_load(self,df):
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
            v["stoich"] = self._str_to_dict(v["stoich"])
            v["element_contribution"] = self._str_to_dict(v["element_contribution"])
        return d
    def load_initiation_subreactions(self):
        create_file = self._modify_initiation_subreactions_for_save(
                self._create_initiation_subreactions())
        df = self._get_manual_curation(
            "initiation_subreactions.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_initiation_subreactions_from_load(df)

    def _create_elongation_subreactions(self):
        return pandas.DataFrame.from_dict(dictionaries.elongation_subreactions.copy()).T.rename_axis('subreaction')
    def _modify_elongation_subreactions_for_save(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
            df.loc[r, "stoich"] = self._dict_to_str(row["stoich"])
        return df
    def _modify_elongation_subreactions_from_load(self,df):
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
            v["stoich"] = self._str_to_dict(v["stoich"])
        return d
    def load_elongation_subreactions(self):
        create_file = self._modify_elongation_subreactions_for_save(
                self._create_elongation_subreactions())
        df = self._get_manual_curation(
            "elongation_subreactions.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_elongation_subreactions_from_load(df)

    def _create_termination_subreactions(self):
        df = pandas.DataFrame.from_dict(dictionaries.termination_subreactions.copy()).T.rename_axis('subreaction')
        df[["element_contribution"]] = df[["element_contribution"]].applymap(
                lambda x: {} if pandas.isnull(x) else x
            )
        return df
    def _modify_termination_subreactions_for_save(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
            df.loc[r, "stoich"] = self._dict_to_str(row["stoich"])
            df.loc[r, "element_contribution"] = self._dict_to_str(
                row["element_contribution"]
            )
        return df
    def _modify_termination_subreactions_from_load(self,df):
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
            v["stoich"] = self._str_to_dict(v["stoich"])
            v["element_contribution"] = self._str_to_dict(v["element_contribution"])
        return d
    def load_termination_subreactions(self):
        create_file = self._modify_termination_subreactions_for_save(
                self._create_termination_subreactions())
        df = self._get_manual_curation(
            "termination_subreactions.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_termination_subreactions_from_load(df)

    def _create_special_trna_subreactions(self):
        df = pandas.DataFrame.from_dict(dictionaries.special_trna_subreactions.copy()).T.rename_axis('subreaction')
        df[["element_contribution"]] = df[["element_contribution"]].applymap(
            lambda x: {} if pandas.isnull(x) else x
        )
        return pandas.DataFrame.from_dict(dictionaries.special_trna_subreactions.copy()).T.rename_axis('subreaction')
    def _modify_special_trna_subreactions_for_save(self,df):

        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
            df.loc[r, "stoich"] = self._dict_to_str(row["stoich"])
            df.loc[r, "element_contribution"] = self._dict_to_str(
                row["element_contribution"]
            )
        return df
    def _modify_special_trna_subreactions_from_load(self,df):
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
            v["stoich"] = self._str_to_dict(v["stoich"])
            v["element_contribution"] = self._str_to_dict(v["element_contribution"])
        return d
    def load_special_trna_subreactions(self):
        create_file = self._modify_special_trna_subreactions_for_save(
                self._create_special_trna_subreactions())
        df = self._get_manual_curation(
            "special_trna_subreactions.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_special_trna_subreactions_from_load(df)

    def _create_excision_machinery(self):
        return pandas.DataFrame.from_dict(dictionaries.excision_machinery.copy()).T.rename_axis('type')
    def _modify_excision_machinery_for_save(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_excision_machinery_from_load(self,df):
        d = df.T.to_dict()
        for k, v in d.items():
            d[k] = {}
            if v['enzymes']:
                d[k]['enzymes'] = v['enzymes'].split(",")
            else:
                d[k]['enzymes'] = []
        return d
    def load_excision_machinery(self):
        create_file = self._modify_excision_machinery_for_save(
                self._create_excision_machinery())
        df = self._get_manual_curation(
            "excision_machinery.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_excision_machinery_from_load(df)

    def _create_special_modifications(self):
        return pandas.DataFrame.from_dict(dictionaries.special_modifications.copy()).T.rename_axis('modification')
    def _modify_special_modifications_for_save(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
            df.loc[r, "stoich"] = self._dict_to_str(row["stoich"])
        return df
    def _modify_special_modifications_from_load(self,df):
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
            v["stoich"] = self._str_to_dict(v["stoich"])
        return d
    def load_special_modifications(self):
        create_file = self._modify_special_modifications_for_save(
                self._create_special_modifications())
        df = self._get_manual_curation(
            "special_modifications.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_special_modifications_from_load(df)

    def _modify_rna_modification_from_load(self,df):
        d = {}
        for mod,row in df.iterrows():
            mods = ['{}_at_{}'.format(mod,i) for i in row['positions'].split(',')]
            for enz in row['enzymes'].split(' AND '):
                if enz not in d: d[enz] = []
                d[enz] = mods
        return d
    def load_rna_modification(self):
        create_file = pandas.DataFrame(columns=[
            'modification','positions','type','enzymes','source'
        ]).set_index('modification')
        df = self._get_manual_curation(
            "rna_modification.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_rna_modification_from_load(df)

#     def _process_rna_modification_targets(self,
#                                df):
#         df = df.reset_index()
#         trna_mod_dict = {}
#         for mod in df.iterrows():
#             mod = mod[1]
#             mod_loc = "%s_at_%s" % (mod["modification"], mod["position"])
#             if mod["bnum"] not in trna_mod_dict:
#                 trna_mod_dict[mod["bnum"]] = {}
#             trna_mod_dict[mod["bnum"]][mod_loc] = 1
#         return trna_mod_dict

    def _create_rna_modification_targets(self):
        return pandas.DataFrame(columns=[
            'bnum',
            'position',
            'modification'
        ]).set_index('bnum')
    def load_rna_modification_targets(self):
        create_file = self._create_rna_modification_targets()
        return self._get_manual_curation(
            "post_transcriptional_modification_of_RNA.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')

    def _create_folding_dict(self):
        return pandas.DataFrame.from_dict(dictionaries.folding_dict.copy()).T.rename_axis('mechanism')
    def _modify_folding_dict_for_save(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_folding_dict_from_load(self,df):
        d = df.T.to_dict()
        for k, v in d.items():
            d[k] = {}
            if v['enzymes']:
                d[k]['enzymes'] = v['enzymes'].split(",")
            else:
                d[k]['enzymes'] = []
        return d
    def load_folding_dict(self):
        create_file = self._modify_folding_dict_for_save(
                self._create_folding_dict())
        df = self._get_manual_curation(
            "folding_dict.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_folding_dict_from_load(df)

    def _create_transcription_subreactions(self):
        return pandas.DataFrame.from_dict(dictionaries.transcription_subreactions.copy()).T.rename_axis('mechanism')
    def _modify_transcription_subreactions_for_save(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
            df.loc[r, "stoich"] = self._dict_to_str(row["stoich"])
        return df
    def _modify_transcription_subreactions_from_load(self,df):
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
            v["stoich"] = self._str_to_dict(v["stoich"])
        return d
    def load_transcription_subreactions(self):
        create_file = self._modify_transcription_subreactions_for_save(
                self._create_transcription_subreactions())
        df = self._get_manual_curation(
            "transcription_subreactions.txt",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_transcription_subreactions_from_load(df)

    def _process_translocation_pathways(self,
                                    df):
        d = {}
        for p in df.index.unique():
            d[p] = {}
            pdf = df.loc[[p]]
            d[p]["keff"] = pdf["keff"][0]
            d[p]["length_dependent_energy"] = pdf["length_dependent_energy"][0]
            d[p]["stoichiometry"] = self._str_to_dict(pdf["stoichiometry"][0])
            d[p]["enzymes"] = {}
            for _, row in pdf.iterrows():
                d[p]["enzymes"][row["enzyme"]] = {
                    "length_dependent": row["length_dependent"],
                    "fixed_keff": row["fixed_keff"],
                }
        return d
    def load_translocation_pathways(self):
        create_file = pandas.DataFrame(columns=
                [
                    "pathway",
                    "keff",
                    "length_dependent_energy",
                    "stoichiometry",
                    "enzyme",
                    "length_dependent",
                    "fixed_keff",
                ]
            ).set_index("pathway")
        return self._process_translocation_pathways(
            self._get_manual_curation(
                "translocation_pathways.txt",
                create_file = create_file,
                no_file_return = create_file,
                sep = '\t'))

    def _modify_lipid_modifications_from_load(self,df):
        return {k:v.split(',') for k,v in df['enzymes'].to_dict().items()}
    def _create_lipid_modifications(self):
        return pandas.DataFrame(columns=[
            'lipid_mod',
            'enzymes'
        ]).set_index('lipid_mod')
    def load_lipid_modifications(self):
        create_file = self._create_lipid_modifications()
        df =  self._get_manual_curation(
                "lipid_modifications.txt",
                create_file = create_file,
                no_file_return = create_file,
                sep = '\t')
        return self._modify_lipid_modifications_from_load(df)

    def _create_subreaction_matrix(self):
        return pandas.DataFrame(columns=[
            'Reaction','Metabolites','Stoichiometry'
        ]).set_index('Reaction')
    def load_subreaction_matrix(self):
        create_file = self._create_subreaction_matrix()
        df =  self._get_manual_curation(
                "subreaction_matrix.txt",
                create_file = create_file,
                no_file_return = create_file,
                sep = '\t')
        return df

    def _create_reaction_matrix(self):
        return pandas.DataFrame(columns=[
            'Reaction', 'Metabolites', 'Stoichiometry'
        ]).set_index('Reaction')
    def load_reaction_matrix(self):
        create_file = self._create_reaction_matrix()
        df =  self._get_manual_curation(
                "reaction_matrix.txt",
                create_file = create_file,
                no_file_return = create_file,
                sep = '\t')
        return df
    
    def _create_stable_RNAs(self):
        return pandas.DataFrame(columns=[
            'id'
        ]).set_index('id')
    def load_stable_RNAs(self):
        create_file = self._create_stable_RNAs()
        df =  self._get_manual_curation(
                "stable_RNAs.txt",
                create_file = create_file,
                no_file_return = create_file,
                sep = '\t')
        return df.index.to_list()

    def _create_orphan_and_spont_reactions(self):
        return pandas.DataFrame(columns=[
            'name', 'description', 'is_reversible', 'is_spontaneous'
        ]).set_index('name')
    def load_orphan_and_spont_reactions(self):
        create_file = self._create_orphan_and_spont_reactions()
        df =  self._get_manual_curation(
                "orphan_and_spont_reactions.txt",
                create_file = create_file,
                no_file_return = create_file,
                sep = '\t')
        return df

    def _create_enz_rxn_assoc_df(self):
        return pandas.DataFrame(columns=[
            'Reaction','Complexes'
        ]).set_index('Reaction')
    def load_enz_rxn_assoc_df(self):
        create_file = self._create_enz_rxn_assoc_df()
        df =  self._get_manual_curation(
                "enzyme_reaction_association.txt",
                create_file = create_file,
                no_file_return = create_file,
                sep = '\t')
        return df

    def _str_to_dict(self,
                    d):
        regex = ":(?=[-]?\d+(?:$|\.))"
        return (
            {re.split(regex, i)[0]: float(re.split(regex, i)[1]) for i in d.split(",")}
            if d
            else {}
        )

    def _dict_to_str(self, d):
        return ",".join(["{}:{}".format(k, v) for k, v in d.items()])



class MECurator(object):
    """MECurator class for integrating additional curation.

    This class integrates additional manual curation that was
    not integrated in Organism.

    Parameters
    ----------
    org : coralme.builder.organism.Organism
        Organism object.
    """

    def __init__(self,
                org):
        self.org = org

#     def curate(self):
#         logging.warning("Integrating manual metabolic reactions")
#         self.modify_metabolic_reactions()
#         logging.warning("Integrating manual complexes")
#         self.add_manual_complexes()

#     def modify_metabolic_reactions(self):
#         m_model = self.org.m_model
#         new_reactions_dict = self.org.reaction_corrections

#         for rxn_id, info in tqdm.tqdm(new_reactions_dict.items(),
#                             'Modifying metabolic reactions with manual curation...',
#                             bar_format = bar_format,
#                             total=len(new_reactions_dict)):
#             if info["reaction"] == "eliminate":
#                 m_model.reactions.get_by_id(rxn_id).remove_from_model()
#             else:
#                 if rxn_id not in m_model.reactions:
#                     rxn = cobra.Reaction(rxn_id)
#                     m_model.add_reaction(rxn)
#                 else:
#                     rxn = m_model.reactions.get_by_id(rxn_id)
#                 if info["reaction"]:
#                     rxn.build_reaction_from_string(info["reaction"])
#                     rxn.name = info["name"]
#                     rxn.gene_reaction_rule = info["gene_reaction_rule"]
#                 if info["gene_reaction_rule"]:
#                     if info["gene_reaction_rule"] == "no_gene":
#                         rxn.gene_reaction_rule = ""
#                     else:
#                         rxn.gene_reaction_rule = info["gene_reaction_rule"]
#                 if info["name"]:
#                     rxn.name = info["name"]

#     def add_manual_complexes(self):
#         manual_complexes = self.org.manual_complexes
#         complexes_df = self.org.complexes_df
#         protein_mod = self.org.protein_mod
#         warn_manual_mod = []
#         warn_replace = []
#         for new_complex, info in tqdm.tqdm(manual_complexes.iterrows(),
#                     'Adding manual curation of complexes...',
#                     bar_format = bar_format,
#                     total=manual_complexes.shape[0]):
#             if info["genes"]:
#                 if new_complex not in complexes_df:
#                     complexes_df = complexes_df.append(
#                         pandas.DataFrame.from_dict(
#                             {new_complex: {"name": "", "genes": "", "source": "Manual"}}
#                         ).T
#                     )
#                 complexes_df.loc[new_complex, "genes"] = info["genes"]
#                 complexes_df.loc[new_complex, "name"] = str(info["name"])
#             if info["mod"]:
#                 mod_complex = (
#                     new_complex
#                     + "".join(
#                         [
#                             "_mod_{}".format(m)
#                             for m in info['mod'].split(' AND ')
#                         ]
#                     )
#                     if info["mod"]
#                     else new_complex
#                 )
#                 if mod_complex in protein_mod.index:
#                     warn_manual_mod.append(mod_complex)
#                     continue
#                 if info["replace"]:
#                     if info["replace"] in protein_mod.index:
#                         protein_mod = protein_mod.drop(info["replace"])
#                     else:
#                         warn_replace.append(mod_complex)
#                 protein_mod = protein_mod.append(
#                     pandas.DataFrame.from_dict(
#                         {
#                             mod_complex: {
#                                 "Core_enzyme": new_complex,
#                                 "Modifications": info["mod"],
#                                 "Source": "Manual",
#                             }
#                         }
#                     ).T
#                 )
#         complexes_df.index.name = "complex"

#         self.org.complexes_df = complexes_df
#         self.org.protein_mod = protein_mod

#         # Warnings
#         if warn_manual_mod or warn_replace:
#             if warn_manual_mod:
#                 self.org.curation_notes['add_manual_complexes'].append({
#                     'msg':'Some modifications in protein_corrections.txt are already in me_builder.org.protein_mod and were skipped.',
#                     'triggered_by':warn_manual_mod,
#                     'importance':'low',
#                     'to_do':'Check whether the protein modification specified in protein_corrections.txt is correct and not duplicated.'})
#             if warn_replace:
#                 self.org.curation_notes['add_manual_complexes'].append({
#                     'msg':'Some modified proteins marked for replacement in protein_corrections.txt are not in me_builder.org.protein_mod. Did nothing.',
#                     'triggered_by':warn_replace,
#                     'importance':'low',
#                     'to_do':'Check whether the marked modified protein in protein_corrections.txt for replacement is correctly defined.'})

    def find_issue_with_query(self,query):
        for k,v in self.org.curation_notes.items():
            coralme.builder.helper_functions.find_issue(query,v)


