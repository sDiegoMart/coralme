import pandas

class MECurator(object):
    
    def __init__(org,
                 config):
        self.org = org
        self.config = config
    
    @property
    def _manual_complexes(self):
        filename = self.directory + "protein_corrections.csv"
        if os.path.isfile(filename):
            return pandas.read_csv(filename, index_col=0).fillna("")
        else:
            self.curation_notes['org._manual_complexes'].append({
                'msg':"No manual complexes file, creating one",
                'importance':'low',
                'to_do':'Fill manual_complexes.csv'})
            df = pandas.DataFrame.from_dict(
                {"complex_id": {}, "name": {}, "genes": {}, "mod": {}, "replace": {}}
            ).set_index("complex_id")
            df.to_csv(filename)
            return df

    @property
    def _sigmas(self):
        filename = self.directory + "sigma_factors.csv"
        if os.path.isfile(filename):
            return pandas.read_csv(
                filename, index_col=0, sep=","
            )
        else:
            return self.get_sigma_factors()

    @property
    def _rpod(self):
        return self.get_rpod()

    @property
    def _m_to_me_mets(self):
        filename = self.directory + "m_to_me_mets.csv"
        if os.path.isfile(filename):
            return pandas.read_csv(filename, index_col=0)
        else:
            self.curation_notes['org._m_to_me_mets'].append({
                'msg':"No m_to_me_mets.csv file found",
                'importance':'low',
                'to_do':'Fill m_to_me_mets.csv'})
            df = pandas.DataFrame.from_dict({"m_name": {}, "me_name": {}})
            df.set_index("m_name")
            df.to_csv(filename)
            return df

    @property
    def _rna_degradosome(self):
        filename = self.directory + "rna_degradosome.csv"
        if os.path.isfile(filename):
            d = (
                pandas.read_csv(filename, index_col=0, sep="\t")
                .fillna("").T
                .to_dict()
            )
            for k, v in d.items():
                d[k] = {}
                if v['enzymes']:
                    d[k]['enzymes'] = v['enzymes'].split(",")
                else:
                    d[k]['enzymes'] = []
        else:
            self.curation_notes['org._rna_degradosome'].append({
                'msg':"No rna_degradosome.csv file found",
                'importance':'low',
                'to_do':'Fill rna_degradosome.csv'})
            d = {'rna_degradosome':{'enzymes':[]}}
            df_save = pandas.DataFrame.from_dict(d).T
            for r, row in df_save.iterrows():
                df_save.loc[r, "enzymes"] = ",".join(row["enzymes"])
            df_save.index.name = "name"
            df_save.to_csv(filename, sep="\t")
        return d

    @property
    def _ribosome_stoich(self):
        filename = self.directory + "ribosomal_proteins.csv"
        if os.path.isfile(filename):
            df = pandas.read_csv(filename, sep="\t", index_col=0).fillna('')
            return self.create_ribosome_stoich(df)
        else:
            self.curation_notes['org._ribosomal_proteins'].append({
                'msg':"No ribosomal_proteins.csv file found",
                'importance':'low',
                'to_do':'Fill ribosomal_proteins.csv'})
            return self.write_ribosome_stoich(filename)

    @property
    def _ribosome_subreactions(self):
        filename = self.directory + "ribosome_subreactions.csv"
        if os.path.isfile(filename):
            d = pandas.read_csv(filename, sep="\t", index_col=0).fillna("").T.to_dict()
            for r, row in d.items():
                d[r]["stoich"] = self.str_to_dict(row["stoich"])
        else:
            self.curation_notes['org._ribosome_subreactions'].append({
                'msg':"No ribosome_subreactions.csv file found",
                'importance':'low',
                'to_do':'Fill ribosome_subreactions.csv'})
            df = pandas.DataFrame.from_dict(dictionaries.ribosome_subreactions).T
            df_save = df.copy()
            for r, row in df_save.iterrows():
                df_save.loc[r, "stoich"] = self.dict_to_str(row["stoich"])
            df_save.to_csv(filename, sep="\t")
            d = df.T.to_dict()
        return d

    @property
    def _generic_dict(self):
        filename = self.directory + "generic_dict.csv"
        if os.path.isfile(filename):
            d = (
                pandas.read_csv(filename, index_col=0, sep="\t")
                .fillna("").T
                .to_dict()
            )
            for k, v in d.items():
                d[k] = {}
                if v['enzymes']:
                    d[k]['enzymes'] = v['enzymes'].split(",")
                else:
                    d[k]['enzymes'] = []
        else:
            self.curation_notes['org._generic_dict'].append({
                'msg':"No generic_dict.csv file found",
                'importance':'low',
                'to_do':'Fill generic_dict.csv'})
            d = dictionaries.generics.copy()
            df_save = pandas.DataFrame.from_dict(d).T
            for r, row in df_save.iterrows():
                df_save.loc[r, "enzymes"] = ",".join(row["enzymes"])
            df_save.index.name = "generic_component"
            df_save.to_csv(filename, sep="\t")
        return d

    @property
    def _rrna_modifications(self):
        filename = self.directory + "rrna_modifications.csv"
        if os.path.isfile(filename):
            d = pandas.read_csv(filename, sep="\t", index_col=0).fillna("").T.to_dict()
            for k, v in d.items():
                v["metabolites"] = self.str_to_dict(v["metabolites"])
        else:
            self.curation_notes['org._rrna_modifications'].append({
                'msg':"No rrna_modifications.csv file found",
                'importance':'low',
                'to_do':'Fill rrna_modifications.csv'})
            df = pandas.DataFrame.from_dict(dictionaries.rrna_modifications).T
            df_save = df.copy()
            for r, row in df_save.iterrows():
                df_save.loc[r, "metabolites"] = self.dict_to_str(row["metabolites"])
            df_save.to_csv(filename, sep="\t")
            d = df.T.to_dict()
        return d

    @property
    def _amino_acid_trna_synthetase(self):
        filename = self.directory + "amino_acid_trna_synthetase.csv"
        if os.path.isfile(filename):
            d = (
                pandas.read_csv(filename, index_col=0, sep="\t")
                .fillna("")["enzyme"]
                .to_dict()
            )
        else:
            self.curation_notes['org._amino_acid_trna_synthetase'].append({
                'msg':"No amino_acid_trna_synthetase.csv file found",
                'importance':'low',
                'to_do':'Fill amino_acid_trna_synthetase.csv'})
            d = dictionaries.amino_acid_trna_synthetase.copy()
            pandas.DataFrame.from_dict({"enzyme": d}).to_csv(filename, sep="\t")
        return d

    @property
    def _peptide_release_factors(self):
        filename = self.directory + "peptide_release_factors.csv"
        if os.path.isfile(filename):
            d = (
                pandas.read_csv(filename, index_col=0, sep="\t")
                .fillna("").T
                .to_dict()
            )
            for k, v in d.items():
                d[k] = {}
                if v['enzyme']:
                    d[k]['enzyme'] = v['enzyme']
                else:
                    d[k]['enzyme'] = ''
        else:
            self.curation_notes['org._peptide_release_factors'].append({
                'msg':"No peptide_release_factors.csv file found",
                'importance':'low',
                'to_do':'Fill peptide_release_factors.csv'})
            d = dictionaries.translation_stop_dict.copy()
            df_save = pandas.DataFrame.from_dict(d).T
            for r, row in df_save.iterrows():
                df_save.loc[r, "enzyme"] = row["enzyme"]
            df_save.index.name = "codon"
            df_save.to_csv(filename, sep="\t")
        return d

    @property
    def _initiation_subreactions(self):
        filename = self.directory + "initiation_subreactions.csv"
        if os.path.isfile(filename):
            d = pandas.read_csv(filename, index_col=0, sep="\t").T.fillna("").to_dict()
            for k, v in d.items():
                v["enzymes"] = v["enzymes"].split(",")
                if "" in v["enzymes"]:
                    v["enzymes"].remove("")
                v["stoich"] = self.str_to_dict(v["stoich"])
                v["element_contribution"] = self.str_to_dict(v["element_contribution"])
        else:
            self.curation_notes['org._initiation_subreactions'].append({
                'msg':"No initiation_subreactions.csv file found",
                'importance':'low',
                'to_do':'Fill initiation_subreactions.csv'})
            df = pandas.DataFrame.from_dict(dictionaries.initiation_subreactions).T
            df_save = df.copy()
            for r, row in df_save.iterrows():
                df_save.loc[r, "enzymes"] = ",".join(row["enzymes"])
                df_save.loc[r, "stoich"] = self.dict_to_str(row["stoich"])
                df_save.loc[r, "element_contribution"] = self.dict_to_str(
                    row["element_contribution"]
                )
            df_save.to_csv(filename, sep="\t")
            d = df.T.to_dict()
        return d

    @property
    def _elongation_subreactions(self):
        filename = self.directory + "elongation_subreactions.csv"
        if os.path.isfile(filename):
            d = pandas.read_csv(filename, index_col=0, sep="\t").T.fillna("").to_dict()
            for k, v in d.items():
                v["enzymes"] = v["enzymes"].split(",")
                if "" in v["enzymes"]:
                    v["enzymes"].remove("")
                v["stoich"] = self.str_to_dict(v["stoich"])
        else:
            self.curation_notes['org._elongation_subreactions'].append({
                'msg':"No elongation_subreactions.csv file found",
                'importance':'low',
                'to_do':'Fill elongation_subreactions.csv'})
            df = pandas.DataFrame.from_dict(dictionaries.elongation_subreactions).T
            df_save = df.copy()
            for r, row in df_save.iterrows():
                df_save.loc[r, "enzymes"] = ",".join(row["enzymes"])
                df_save.loc[r, "stoich"] = self.dict_to_str(row["stoich"])
            df_save.to_csv(filename, sep="\t")
            d = df.T.to_dict()
        return d

    @property
    def _termination_subreactions(self):
        filename = self.directory + "termination_subreactions.csv"
        if os.path.isfile(filename):
            d = pandas.read_csv(filename, index_col=0, sep="\t").T.fillna("").to_dict()
            for k, v in d.items():
                v["enzymes"] = v["enzymes"].split(",")
                if "" in v["enzymes"]:
                    v["enzymes"].remove("")
                v["stoich"] = self.str_to_dict(v["stoich"])
                v["element_contribution"] = self.str_to_dict(v["element_contribution"])
        else:
            self.curation_notes['org._termination_subreactions'].append({
                'msg':"No termination_subreactions.csv file found",
                'importance':'low',
                'to_do':'Fill termination_subreactions.csv'})
            df = pandas.DataFrame.from_dict(dictionaries.termination_subreactions).T
            df[["element_contribution"]] = df[["element_contribution"]].applymap(
                lambda x: {} if pandas.isnull(x) else x
            )
            df_save = df.copy()
            for r, row in df_save.iterrows():
                df_save.loc[r, "enzymes"] = ",".join(row["enzymes"])
                df_save.loc[r, "stoich"] = self.dict_to_str(row["stoich"])
                df_save.loc[r, "element_contribution"] = self.dict_to_str(
                    row["element_contribution"]
                )
            df_save.to_csv(filename, sep="\t")
            d = df.T.to_dict()
        return d

    @property
    def _special_trna_subreactions(self):
        filename = self.directory + "special_trna_subreactions.csv"
        if os.path.isfile(filename):
            d = pandas.read_csv(filename, index_col=0, sep="\t").T.fillna("").to_dict()
            for k, v in d.items():
                v["enzymes"] = v["enzymes"].split(",")
                if "" in v["enzymes"]:
                    v["enzymes"].remove("")
                v["stoich"] = self.str_to_dict(v["stoich"])
                v["element_contribution"] = self.str_to_dict(v["element_contribution"])
        else:
            self.curation_notes['org._special_trna_subreactions'].append({
                'msg':"No special_trna_subreactions.csv file found",
                'importance':'low',
                'to_do':'Fill special_trna_subreactions.csv'})
            df = pandas.DataFrame.from_dict(dictionaries.special_trna_subreactions).T
            df_save = df.copy()
            for r, row in df_save.iterrows():
                df_save.loc[r, "enzymes"] = ",".join(row["enzymes"])
                df_save.loc[r, "stoich"] = self.dict_to_str(row["stoich"])
                df_save.loc[r, "element_contribution"] = self.dict_to_str(
                    row["element_contribution"]
                )
            df_save.to_csv(filename, sep="\t")
            d = df.T.to_dict()
        return d

    @property
    def _excision_machinery(self):
        filename = self.directory + "excision_machinery.csv"
        if os.path.isfile(filename):
            d = (
                pandas.read_csv(filename, index_col=0, sep="\t")
                .fillna("").T
                .to_dict()
            )
            for k, v in d.items():
                d[k] = {}
                if v['enzymes']:
                    d[k]['enzymes'] = v['enzymes'].split(",")
                else:
                    d[k]['enzymes'] = []
        else:
            self.curation_notes['org._excision_machinery'].append({
                'msg':"No excision_machinery.csv file found",
                'importance':'low',
                'to_do':'Fill excision_machinery.csv'})
            d = dictionaries.excision_machinery.copy()
            df_save = pandas.DataFrame.from_dict(d).T
            for r, row in df_save.iterrows():
                df_save.loc[r, "enzymes"] = ",".join(row["enzymes"])
            df_save.index.name = "mechanism"
            df_save.to_csv(filename, sep="\t")
        return d

    @property
    def _special_modifications(self):
        filename = self.directory + "special_modifications.csv"
        if os.path.isfile(filename):
            d = pandas.read_csv(filename, index_col=0, sep="\t").T.fillna("").to_dict()
            for k, v in d.items():
                v["enzymes"] = v["enzymes"].split(",")
                if "" in v["enzymes"]:
                    v["enzymes"].remove("")
                v["stoich"] = self.str_to_dict(v["stoich"])
        else:
            self.curation_notes['org._special_modifications'].append({
                'msg':"No special_modifications.csv file found",
                'importance':'low',
                'to_do':'Fill special_modifications.csv'})
            df = pandas.DataFrame.from_dict(dictionaries.special_modifications).T
            df_save = df.copy()
            for r, row in df_save.iterrows():
                df_save.loc[r, "enzymes"] = ",".join(row["enzymes"])
                df_save.loc[r, "stoich"] = self.dict_to_str(row["stoich"])
            df_save.to_csv(filename, sep="\t")
            d = df.T.to_dict()
        return d

    @property
    def _trna_modification(self):
        filename = self.directory + "trna_modification.csv"
        if os.path.isfile(filename):
            d = pandas.read_csv(filename, index_col=0, sep="\t").T.fillna("").to_dict()
            for k, v in d.items():
                v["enzymes"] = v["enzymes"].split(",")
                if "" in v["enzymes"]:
                    v["enzymes"].remove("")
                v["stoich"] = self.str_to_dict(v["stoich"])
                v["carriers"] = self.str_to_dict(v["carriers"])
        else:
            self.curation_notes['org._trna_modification'].append({
                'msg':"No trna_modification.csv file found",
                'importance':'low',
                'to_do':'Fill trna_modification.csv'})
            df = pandas.DataFrame.from_dict(dictionaries.trna_modification).T
            df[["carriers"]] = df[["carriers"]].applymap(
                lambda x: {} if pandas.isnull(x) else x
            )
            df_save = df.copy()
            for r, row in df_save.iterrows():
                df_save.loc[r, "enzymes"] = ",".join(row["enzymes"])
                df_save.loc[r, "stoich"] = self.dict_to_str(row["stoich"])
                df_save.loc[r, "carriers"] = self.dict_to_str(row["carriers"])
            df_save.to_csv(filename, sep="\t")
            d = df.T.to_dict()
        return d

    @property
    def _trna_modification_targets(self):
        filename = self.directory + "post_transcriptional_modification_of_tRNA.csv"
        if os.path.isfile(filename):
            df = pandas.read_csv(filename, delimiter="\t")
        else:
            self.curation_notes['org._trna_modification_targets'].append({
                'msg':"No post_transcriptional_modification_of_tRNA.csv file found",
                'importance':'low',
                'to_do':'Fill post_transcriptional_modification_of_tRNA.csv'})
            df = pandas.DataFrame.from_dict(
                {"bnum": {}, "position": {}, "modification": {}})

            df.to_csv(filename, sep="\t")
        trna_mod_dict = {}
        for mod in df.iterrows():
            mod = mod[1]
            mod_loc = "%s_at_%s" % (mod["modification"], mod["position"])
            if mod["bnum"] not in trna_mod_dict:
                trna_mod_dict[mod["bnum"]] = {}
            trna_mod_dict[mod["bnum"]][mod_loc] = 1
        return trna_mod_dict

    @property
    def _folding_dict(self):
        filename = self.directory + "folding_dict.csv"
        if os.path.isfile(filename):
            d = (
                pandas.read_csv(filename, index_col=0, sep="\t")
                .fillna("").T
                .to_dict()
            )
            for k, v in d.items():
                d[k] = {}
                if v['enzymes']:
                    d[k]['enzymes'] = v['enzymes'].split(",")
                else:
                    d[k]['enzymes'] = []
        else:
            self.curation_notes['org._folding_dict'].append({
                'msg':"No folding_dict.csv file found",
                'importance':'low',
                'to_do':'Fill folding_dict.csv'})
            d = dictionaries.folding_dict.copy()
            df_save = pandas.DataFrame.from_dict(d).T
            for r, row in df_save.iterrows():
                df_save.loc[r, "enzymes"] = ",".join(row["enzymes"])
            df_save.index.name = "mechanism"
            df_save.to_csv(filename, sep="\t")
        return d

    @property
    def _transcription_subreactions(self):
        filename = self.directory + "transcription_subreactions.csv"
        if os.path.isfile(filename):
            d = pandas.read_csv(filename, index_col=0, sep="\t").T.fillna("").to_dict()
            for k, v in d.items():
                v["enzymes"] = v["enzymes"].split(",")
                if "" in v["enzymes"]:
                    v["enzymes"].remove("")
                v["stoich"] = self.str_to_dict(v["stoich"])
        else:
            self.curation_notes['org._transcription_subreactions'].append({
                'msg':"No transcription_subreactions.csv file found",
                'importance':'low',
                'to_do':'Fill transcription_subreactions.csv'})
            df = pandas.DataFrame.from_dict(dictionaries.transcription_subreactions).T
            df_save = df.copy()
            for r, row in df_save.iterrows():
                df_save.loc[r, "enzymes"] = ",".join(row["enzymes"])
                df_save.loc[r, "stoich"] = self.dict_to_str(row["stoich"])
            df_save.to_csv(filename, sep="\t")
            d = df.T.to_dict()
        return d

    @property
    def _translocation_pathways(self):
        filename = self.directory + "translocation_pathways.csv"
        if os.path.isfile(filename):
            df = pandas.read_csv(filename, index_col=0, delimiter="\t").fillna("")
            return self.read_translocation_pathways(df)
        else:
            self.curation_notes['org._translocation_pathways'].append({
                'msg':"No translocation_pathways.csv file found",
                'importance':'low',
                'to_do':'Fill translocation_pathways.csv'})
            columns = [
                "keff",
                "length_dependent_energy",
                "stoichiometry",
                "enzyme",
                "length_dependent",
                "fixed_keff",
            ]
            df = pandas.DataFrame(columns=columns)
            df.index.name = 'pathway'
            df.to_csv(filename, sep="\t")
            return dict()



    def _get_manual_curation(self,
                             filename,
                             create_file=None,
                             no_file_return=pandas.DataFrame(),
                             sep = '\t'):
        filepath = self.directory + filename
        if os.path.isfile(filepath):
            return pandas.read_csv(filepath, index_col=0,sep=sep)
        
        if create_file is not None:
            create_file.to_csv(filepath)
            
        self.curation_notes['org._get_manual_curation'].append({
            'msg':'No {} file found'.format(filename),
            'importance':'low',
            'to_do':'Fill in {}'.format(filepath)
        })
        return no_file_return
    
    ## TODO: implement this
    def load_protein_location(self):
        return self._get_manual_curation(
             "peptide_compartment_and_pathways.csv").to_dict() 
        filename = self.directory + "peptide_compartment_and_pathways.csv"
        if os.path.isfile(filename):
            return pandas.read_csv(filename, index_col = 0, delimiter = "\t")
        else:
            self.curation_notes['org._protein_location'].append({
                'msg':"No peptide_compartment_and_pathways.csv file found",
                'importance':'low',
                'to_do':'Fill peptide_compartment_and_pathways.csv'})
            columns = [
                'Complex',
                'Complex_compartment',
                'Protein',
                'Protein_compartment',
                'translocase_pathway',
                ]
            pandas.DataFrame(columns = columns).set_index('Complex').to_csv(filename, sep = "\t")
            return self.get_protein_location()

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

    def load_subsystem_classification(self):
        if self.is_reference:
            return None

        subsystems = set(r.subsystem for r in self.m_model.reactions if r.subsystem)
        create_file = {}
        for s in subsystems:
            create_file[s] = {}
            create_file[s]["central_CE"] = 0
            create_file[s]["central_AFN"] = 0
            create_file[s]["intermediate"] = 0
            create_file[s]["secondary"] = 0
            create_file[s]["other"] = 1
        create_file = pandas.DataFrame.from_dict(create_file).T

        d = {}
        df = self._get_manual_curation(
            "subsystem_classification.csv",
            create_file = create_file,
            no_file_return = pandas.DataFrame())

        for c in df.columns:
            for s in df[df[c] == 1].index:
                d[s] = c
        return d