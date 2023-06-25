import os
import re

import logging
log = logging.getLogger(__name__)

import tqdm
bar_format = '{desc:<75}: {percentage:.1f}%|{bar}| {n_fmt:>5}/{total_fmt:>5} [{elapsed}<{remaining}]'

import coralme
import cobra
import pandas
import json
import copy

cur_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(cur_dir, 'column_format.json'), 'r') as f:
    column_format = json.load(f)

class CurationInfo(object):
    def __init__(self,
                 id,
                 org,
                 config = {},
                 file = ''):
        self.id = id
        if file:
            self.file = file
        else:
            self.file = id + '.txt'
        self.directory = org.directory
        self.org = org
        self.config = config
        self.data = self.load()
        self.org.__setattr__(id,copy.copy(self.data))
    
    def read(self):
        if self.config["pathtype"] == 'absolute':
            self.filepath = self.file
        else:
            self.filepath = self.directory + self.file
        if os.path.isfile(self.filepath):
            return pandas.read_csv(self.filepath,
                                   index_col=0,
                                   sep=self.config["sep"],
                                   comment='#',
                                   skip_blank_lines=True).fillna("")
        return None
    
    def load(self):
        self.data = self.read()
        if self.data is None:
            self.org.curation_notes['org._get_manual_curation'].append({
                'msg':'No {} file found'.format(self.id),
                'importance':'low',
                'to_do':'Fill in {}'.format(self.filepath)
            })
            self.data = self._modify_for_create(self.config["create_file"])
            self.data.to_csv(self.filepath,sep=self.config["sep"])
        self.data = self._modify_from_load()
        return self.data
    
    def save(self):
        return
    def _modify_from_load(self):
        return self.data
    def _modify_for_save(self):
        return self.data
    def _modify_for_create(self,df):
        return df
    
    @property
    def columns(self):
        return column_format[self.file]

class ReactionCorrections(CurationInfo):
    def __init__(self,
                 org,
                 id = "reaction_corrections",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : ',',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_from_load(self):
        return self.data.T.to_dict()
    
class ProteinLocation(CurationInfo):
    def __init__(self,
                 org,
                 id = "protein_location",
                 config={},
                 file="peptide_compartment_and_pathways.txt"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    
class TranslocationMultipliers(CurationInfo):
    def __init__(self,
                 org,
                 id = "translocation_multipliers",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        create_file = pandas.DataFrame()
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : ',',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_from_load(self):
        return self.data.to_dict()
    
class LipoproteinPrecursors(CurationInfo):
    def __init__(self,
                 org,
                 id = "lipoprotein_precursors",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : ',',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_from_load(self):
        return self.data.to_dict()[self.columns[1]]
    
class CleavedMethionine(CurationInfo):
    def __init__(self,
                 org,
                 id = "cleaved_methionine",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_from_load(self):
        return self.data.index.to_list()
    
class ManualComplexes(CurationInfo):
    def __init__(self,
                 org,
                 id = "manual_complexes",
                 config={},
                 file="protein_corrections.txt"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : ',',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    
class Sigmas(CurationInfo):
    def __init__(self,
                 org,
                 id = "sigmas",
                 config={},
                 file="sigma_factors.txt"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : ',',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)

class RhoIndependent(CurationInfo):
    def __init__(self,
                 org,
                 id = "rho_independent",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_from_load(self):
        return self.data.index.to_list()
    
class RNADegradosome(CurationInfo):
    def __init__(self,
                 org,
                 id = "rna_degradosome",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    
    def _modify_from_load(self):
        return {"rna_degradosome" : {"enzymes" : self.data.index.to_list()}}
    
class RNAModificationMachinery(CurationInfo):
    def __init__(self,
                 org,
                 id = "rna_modification_df",
                 config={},
                 file="rna_modification.txt"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    
    def _modify_from_load(self):
        return self.data.astype(str)
    
class RNAModificationTargets(CurationInfo):
    def __init__(self,
                 org,
                 id = "rna_modification_targets",
                 config={},
                 file="post_transcriptional_modification_of_RNA.txt"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    
class EnzymeReactionAssociation(CurationInfo):
    def __init__(self,
                 org,
                 id = "enz_rxn_assoc_df",
                 config={},
                 file="enzyme_reaction_association.txt"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    
class MEMetabolites(CurationInfo):
    def __init__(self,
                 org,
                 id = "me_mets",
                 config={},
                 file="me_metabolites.txt"):
        if not file:
            file = id + ".txt"        
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    
class SubreactionMatrix(CurationInfo):
    def __init__(self,
                 org,
                 id = "subreaction_matrix",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    
class ReactionMatrix(CurationInfo):
    def __init__(self,
                 org,
                 id = "reaction_matrix",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    
class OrphanSpontReactions(CurationInfo):
    def __init__(self,
                 org,
                 id = "orphan_and_spont_reactions",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    
class SubsystemClassification(CurationInfo):
    def __init__(self,
                 org,
                 id = "subsystem_classification",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_from_load(self):
        subsystems = set(r.subsystem for r in self.org.m_model.reactions if r.subsystem)
        df = self.data
        d = {}
        for c in df.columns:
            for s in df[df[c] == 1].index:
                d[s] = c
        return d
    
class TranslocationPathways(CurationInfo):
    def __init__(self,
                 org,
                 id = "translocation_pathways",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = {}
        for p in df.index.unique():
            d[p] = {}
            pdf = df.loc[[p]]
            d[p]["enzymes"] = {}
            for _, row in pdf.iterrows():
                d[p]["enzymes"][row["enzyme"]] = {
                    "length_dependent": row["length_dependent"],
                    "fixed_keff": row["fixed_keff"],
                }
        return d
    
class LipodModifications(CurationInfo):
    def __init__(self,
                 org,
                 id = "lipid_modifications",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_from_load(self):
        df = self.data
        return {k:v.split(',') for k,v in df['enzymes'].to_dict().items()}
    
class StableRNAs(CurationInfo):
    def __init__(self,
                 org,
                 id = "stable_RNAs",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame(columns = self.columns).set_index(self.columns[0])
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_from_load(self):
        return self.data.index.to_list()

        
class RibosomeStoich(CurationInfo):
    def __init__(self,
                 org,
                 id = "ribosome_stoich",
                 config={},
                 file="ribosomal_proteins.txt"):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            {"proteins": {"30S": "generic_16s_rRNAs",
                          "50S": "generic_5s_rRNAs,generic_23s_rRNAs"}}).rename_axis('subunit')
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_from_load(self):
        df = self.data
        from copy import deepcopy
        ribosome_stoich = deepcopy(coralme.builder.dictionaries.ribosome_stoich)
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
    
class RibosomeSubreactions(CurationInfo):
    def __init__(self,
                 org,
                 id = "ribosome_subreactions",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(coralme.builder.dictionaries.ribosome_subreactions.copy()).T.rename_axis('subreaction')\
                .drop("stoich",axis=1).drop("num_mods",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)    
    def _modify_from_load(self):
        return self.data.T.to_dict()
        
class GenericDict(CurationInfo):
    def __init__(self,
                 org,
                 id = "generic_dict",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(coralme.builder.dictionaries.generics.copy()).T.rename_axis('generic_component')
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            d[k] = {}
            if v['enzymes']:
                d[k]['enzymes'] = v['enzymes'].split(",")
            else:
                d[k]['enzymes'] = []
        return d
    
class AminoacidtRNASynthetase(CurationInfo):
    def __init__(self,
                 org,
                 id = "amino_acid_trna_synthetase",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            {"enzyme": coralme.builder.dictionaries.amino_acid_trna_synthetase.copy()}).rename_axis('amino_acid')
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_from_load(self):
        return self.data.to_dict()['enzyme']
    def _modify_for_create(self,df):
        return df
    
class PeptideReleaseFactors(CurationInfo):
    def __init__(self,
                 org,
                 id = "peptide_release_factors",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.translation_stop_dict.copy()).T.rename_axis('release_factor')
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzyme"] = row["enzyme"]
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            d[k] = {}
            if v['enzyme']:
                d[k]['enzyme'] = v['enzyme']
            else:
                d[k]['enzyme'] = ''
        return d
    
class InitiationSubreactions(CurationInfo):
    def __init__(self,
                 org,
                 id = "initiation_subreactions",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.initiation_subreactions.copy()).T.rename_axis('subreaction')\
            .drop("stoich",axis=1).drop("element_contribution",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
        return d
    
class ElongationSubreactions(CurationInfo):
    def __init__(self,
                 org,
                 id = "elongation_subreactions",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.elongation_subreactions.copy()).T.rename_axis('subreaction')\
            .drop("stoich",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
        return d
    
class TerminationSubreactions(CurationInfo):
    def __init__(self,
                 org,
                 id = "termination_subreactions",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.termination_subreactions.copy()).T.rename_axis('subreaction')\
            .drop("stoich",axis=1).drop("element_contribution",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
        return d
    
class SpecialtRNASubreactions(CurationInfo):
    def __init__(self,
                 org,
                 id = "special_trna_subreactions",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.special_trna_subreactions.copy()).T.rename_axis('subreaction')\
            .drop("stoich",axis=1).drop("element_contribution",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
        return d
    
class TranscriptionSubreactions(CurationInfo):
    def __init__(self,
                 org,
                 id = "transcription_subreactions",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.transcription_subreactions.copy()).T.rename_axis('mechanism')\
            .drop("stoich",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
        return d
    
class SpecialModifications(CurationInfo):
    def __init__(self,
                 org,
                 id = "special_modifications",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(coralme.builder.dictionaries.special_modifications.copy()).T.rename_axis('modification')\
            .drop("stoich",axis=1)
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            v["enzymes"] = v["enzymes"].split(",")
            if "" in v["enzymes"]:
                v["enzymes"].remove("")
        return d
    
class ExcisionMachinery(CurationInfo):
    def __init__(self,
                 org,
                 id = "excision_machinery",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.excision_machinery.copy()).T.rename_axis('type')
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            d[k] = {}
            if v['enzymes']:
                d[k]['enzymes'] = v['enzymes'].split(",")
            else:
                d[k]['enzymes'] = []
        return d

class FoldingDict(CurationInfo):
    def __init__(self,
                 org,
                 id = "folding_dict",
                 config={},
                 file=""):
        if not file:
            file = id + ".txt"
        self.file = file
        create_file = pandas.DataFrame.from_dict(
            coralme.builder.dictionaries.folding_dict.copy()).T.rename_axis('mechanism')
        if not config:
            config = {
                     "create_file" : create_file,
                     "sep" : '\t',
                     "pathtype" : 'relative'
                }
        super().__init__(id,
                        org,
                        config = config,
                        file = file)
    def _modify_for_create(self,df):
        df = df.copy()
        for r, row in df.iterrows():
            df.loc[r, "enzymes"] = ",".join(row["enzymes"])
        return df
    def _modify_from_load(self):
        df = self.data
        d = df.T.to_dict()
        for k, v in d.items():
            d[k] = {}
            if v['enzymes']:
                d[k]['enzymes'] = v['enzymes'].split(",")
            else:
                d[k]['enzymes'] = []
        return d
        
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
        self.org.manual_curation = cobra.core.dictlist.DictList()
        self.org.manual_curation.append(ReactionCorrections(self.org))
        self.org.manual_curation.append(ProteinLocation(self.org))
        self.org.manual_curation.append(TranslocationMultipliers(self.org))
        self.org.manual_curation.append(LipoproteinPrecursors(self.org))
        self.org.manual_curation.append(CleavedMethionine(self.org))
        self.org.manual_curation.append(ManualComplexes(self.org))
        self.org.manual_curation.append(Sigmas(self.org))
        self.org.manual_curation.append(RhoIndependent(self.org))
        self.org.manual_curation.append(RNADegradosome(self.org))
        self.org.manual_curation.append(RNAModificationMachinery(self.org))
        self.org.manual_curation.append(RNAModificationTargets(self.org))
        self.org.manual_curation.append(EnzymeReactionAssociation(self.org))
        self.org.manual_curation.append(MEMetabolites(self.org))
        self.org.manual_curation.append(SubreactionMatrix(self.org))
        self.org.manual_curation.append(ReactionMatrix(self.org))
        self.org.manual_curation.append(OrphanSpontReactions(self.org))
        self.org.manual_curation.append(SubsystemClassification(self.org))
        self.org.manual_curation.append(TranslocationPathways(self.org))
        self.org.manual_curation.append(LipodModifications(self.org))
        self.org.manual_curation.append(StableRNAs(self.org))
        self.org.manual_curation.append(RibosomeStoich(self.org))
        self.org.manual_curation.append(RibosomeSubreactions(self.org))
        self.org.manual_curation.append(GenericDict(self.org))
        self.org.manual_curation.append(AminoacidtRNASynthetase(self.org))
        self.org.manual_curation.append(PeptideReleaseFactors(self.org))
        self.org.manual_curation.append(InitiationSubreactions(self.org))
        self.org.manual_curation.append(ElongationSubreactions(self.org))
        self.org.manual_curation.append(TerminationSubreactions(self.org))
        self.org.manual_curation.append(TranscriptionSubreactions(self.org))
        self.org.manual_curation.append(SpecialtRNASubreactions(self.org))
        self.org.manual_curation.append(SpecialModifications(self.org))
        self.org.manual_curation.append(ExcisionMachinery(self.org))
        self.org.manual_curation.append(FoldingDict(self.org))


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

    def find_issue_with_query(self,query):
        for k,v in self.org.curation_notes.items():
            coralme.builder.helper_functions.find_issue(query,v)

            
def _str_to_dict(d):
    regex = ":(?=[-]?\d+(?:$|\.))"
    return (
        {re.split(regex, i)[0]: float(re.split(regex, i)[1]) for i in d.split(",")}
        if d
        else {}
    )

def _dict_to_str(d):
    return ",".join(["{}:{}".format(k, v) for k, v in d.items()])

