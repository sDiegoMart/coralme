import pandas
import os
from coralme.builder import dictionaries
import re

class MECurator(object):
    
    def __init__(self,
                 builder):
        self.builder = builder
        self.directory = builder.org.directory
        self.curation_notes = builder.org.curation_notes
        self.is_reference = builder.org.is_reference

    def _get_manual_curation(self,
                             filename,
                             create_file=None,
                             no_file_return=pandas.DataFrame(),
                             sep = '\t'):
        filepath = self.directory + filename
        if os.path.isfile(filepath):
            return pandas.read_csv(filepath, index_col=0,sep=sep).fillna("")
        
        if create_file is not None:
            create_file.to_csv(filepath,sep=sep)
            
        self.curation_notes['org._get_manual_curation'].append({
            'msg':'No {} file found'.format(filename),
            'importance':'low',
            'to_do':'Fill in {}'.format(filepath)
        })
        return no_file_return
    
    def load_protein_location(self):
        df = pandas.DataFrame(columns = [
                'Complex',
                'Complex_compartment',
                'Protein',
                'Protein_compartment',
                'translocase_pathway',
                ]).set_index('Complex')
        return self._get_manual_curation(
             "peptide_compartment_and_pathways.csv",
             create_file = df,
             no_file_return = df)
        filename = self.directory + "peptide_compartment_and_pathways.csv"

    def load_translocation_multipliers(self):
        return self._get_manual_curation(
             "translocation_multipliers.csv").to_dict()

    def load_lipoprotein_precursors(self):
        return self._get_manual_curation(
            "lipoprotein_precursors.csv",
            no_file_return = pandas.DataFrame(columns=['gene'])).to_dict()["gene"]

    def load_cleaved_methionine(self):
        return self._get_manual_curation(
            "cleaved_methionine.csv",
            create_file = pandas.DataFrame.from_dict({'cleaved_methionine_genes':{}}).set_index('cleaved_methionine_genes'),
            no_file_return = list())

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

        subsystems = set(r.subsystem for r in self.builder.org.m_model.reactions if r.subsystem)
        create_file = self._create_subsystem_classification(subsystems)

        df = self._get_manual_curation(
            "subsystem_classification.csv",
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
            "protein_corrections.csv",
            create_file = create_file,
            no_file_return = create_file,
            sep = ',')

    def load_sigmas(self):
        create_file = pandas.DataFrame(columns = [
                'sigma','complex','genes','name'
            ]).set_index('sigma')
        return self._get_manual_curation(
            "sigma_factors.csv",
            create_file = create_file,
            no_file_return = create_file,
            sep = ',')

    def load_m_to_me_mets(self):
        create_file = pandas.DataFrame(columns = [
                'm_name','me_name'
            ]).set_index('m_name')
        return self._get_manual_curation(
            "m_to_me_mets.csv",
            create_file = create_file,
            no_file_return = create_file,
            sep = ',')
    
    def load_rna_degradosome(self):
        create_file = pandas.DataFrame(columns = [
                'rna_degradosome'
            ]).set_index('rna_degradosome')
        df = self._get_manual_curation(
            "rna_degradosome.csv",
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
            "ribosomal_proteins.csv",
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
            "ribosome_subreactions.csv",
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
            "generic_dict.csv",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_generic_dict_from_load(df)
    
    def _create_rrna_modifications(self):
        return pandas.DataFrame.from_dict(dictionaries.rrna_modifications.copy()).T.rename_axis('modification')
    def _modify_rrna_modifications_for_save(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "metabolites"] = self._dict_to_str(row["metabolites"])
        return df
    def _modify_rrna_modifications_from_load(self,df):
        d = df.T.to_dict()
        for k, v in d.items():
            v["metabolites"] = self._str_to_dict(v["metabolites"])
        return d
    def load_rrna_modifications(self):
        create_file = self._modify_rrna_modifications_for_save(
                self._create_rrna_modifications()) 
        df = self._get_manual_curation(
            "rrna_modifications.csv",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_rrna_modifications_from_load(df)
    
    def _create_amino_acid_trna_synthetase(self):
        return pandas.DataFrame.from_dict(
            {"enzyme": dictionaries.amino_acid_trna_synthetase.copy()}).rename_axis('amino_acid')
    def load_amino_acid_trna_synthetase(self):
        create_file = self._create_amino_acid_trna_synthetase()
        return self._get_manual_curation(
                "amino_acid_trna_synthetase.csv",
                create_file = create_file,
                no_file_return = create_file,
                sep = '\t').to_dict()
    
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
            "peptide_release_factors.csv",
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
            "initiation_subreactions.csv",
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
            "elongation_subreactions.csv",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_elongation_subreactions_from_load(df)
    
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
            "special_trna_subreactions.csv",
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
            "excision_machinery.csv",
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
            "special_modifications.csv",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_special_modifications_from_load(df)
    
    def _create_trna_modification(self):
        df = pandas.DataFrame.from_dict(dictionaries.trna_modification.copy()).T.rename_axis('modification')
        df[["carriers"]] = df[["carriers"]].applymap(
                lambda x: {} if pandas.isnull(x) else x
            )
        return df
    def _modify_trna_modification_for_save(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
            df.loc[r, "stoich"] = self._dict_to_str(row["stoich"])
            df.loc[r, "carriers"] = self._dict_to_str(row["carriers"])
        return df
    def _modify_trna_modification_from_load(self,df):
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
            v["stoich"] = self._str_to_dict(v["stoich"])
            v["carriers"] = self._str_to_dict(v["carriers"])
        return d
    def load_trna_modification(self):
        create_file = self._modify_trna_modification_for_save(
                self._create_trna_modification()) 
        df = self._get_manual_curation(
            "trna_modification.csv",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._modify_trna_modification_from_load(df)

    def _process_trna_modification(self,
                               df):
        trna_mod_dict = {}
        for mod in df.iterrows():
            mod = mod[1]
            mod_loc = "%s_at_%s" % (mod["modification"], mod["position"])
            if mod["bnum"] not in trna_mod_dict:
                trna_mod_dict[mod["bnum"]] = {}
            trna_mod_dict[mod["bnum"]][mod_loc] = 1
        return trna_mod_dict
    
    def _create_trna_modification(self):
        return pandas.DataFrame(columns=[
            'bnum',
            'position',
            'modification'
        ]).set_index('bnum')
    
    def load_trna_modification_targets(self):
        create_file = self._create_trna_modification()
        df = self._get_manual_curation(
            "post_transcriptional_modification_of_tRNA.csv",
            create_file = create_file,
            no_file_return = create_file,
            sep = '\t')
        return self._process_trna_modification(df)

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
            "folding_dict.csv",
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
            "transcription_subreactions.csv",
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
            d[p]["stoichiometry"] = self.str_to_dict(pdf["stoichiometry"][0])
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
                "translocation_pathways.csv",
                create_file = create_file,
                no_file_return = create_file,
                sep = '\t'))
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