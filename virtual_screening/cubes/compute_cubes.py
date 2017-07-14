
from __future__ import print_function
import sys, os
import json,requests
from requests_toolbelt import MultipartEncoder
import urllib.parse as parse

import tempfile

import time
import random

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from openeye import oechem
from openeye import oegraphsim
from openeye.oegraphsim import OEFPAtomType_AtomicNumber as AtmNum, OEFPAtomType_Aromaticity as Arom, OEFPAtomType_Chiral as Chiral, OEFPAtomType_FormalCharge as FCharge, OEFPAtomType_HvyDegree as HvyDeg, OEFPAtomType_Hybridization as Hyb, OEFPAtomType_EqHalogen as EqHalo, OEFPAtomType_HCount as HCount 
from openeye.oegraphsim import OEFPBondType_BondOrder as Order, OEFPBondType_Chiral as ChiralB
from openeye import oeshape
from openeye import oeomega
from openeye import oemolprop
from openeye import oefastrocs 

try:
    from xmlrpclib import ServerProxy, Binary, Fault
except ImportError: # python 3
    from xmlrpc.client import ServerProxy, Binary, Fault

from cubes.vs_classes import VirtualScreeningData, ActiveList, ObjectInputPort, ObjectOutputPort

#from floe.api.parameter import (IntegerParameter, DataSetInputParameter, FileOutputParameter, FileInputParameter,
#                                DataSetOutputParameter, BaseParameter, ParameterGroup,
#                                DecimalParameter, StringParameter)

from floe.api.ports import (InputPort, OutputPort, Port, MoleculeInputPort,
                            MoleculeOutputPort, BinaryMoleculeInputPort, BinaryOutputPort,
                            MoleculeSerializerMixin, BinaryInputPort)

from floe.api import (parameter, ParallelOEMolComputeCube, OEMolComputeCube, SourceCube,
                      ComputeCube, ParallelComputeCube,
                      MoleculeOutputPort)


class AccuMolList(ComputeCube):
    """
    A compute Cube that gets Molecules and assemble them in a list
    """

    classification = [["Compute", "Concat"]]

    intake = MoleculeInputPort('intake')
    success = ObjectOutputPort('success')

    def begin(self):
        self.mol_list = list()

    def process(self, mol, port):
        self.mol_list.append(mol)

    def end(self):
        self.success.emit(self.mol_list)
        
class IndexGenerator(ComputeCube):
    """

    """

    classification = [["Compute", "Index"]]

    intake = ObjectInputPort('intake')
    success = ObjectOutputPort('success')
    tags =[['IndexGenerator_tags']]
    title = 'IndexGenerator_Title' 

    def begin(self):
        self.log.info("Generator Cube running")

    def process(self, data, port):
        self.act_list = data
        pass

    def end(self):
        total = len(self.act_list)
        nb_index = total // 3
        iteration = 10
        
        for idx in range (iteration):
            self.baitset = list()
            i = 0

            index_set = set()
            while i < nb_index :
                index = random.randint(0, total - 1)
                if index in index_set:
                    continue
                else : 
                    index_set.add(index)
                    self.baitset.append(index)
                    i += 1

            self.baitset.sort()
            self.success.emit((idx, self.baitset))

class PrepareRanking(ComputeCube):

    url = parameter.StringParameter('url', default="", help_text="Url of the Restful FastROCS Server for the request")

    act_input = ObjectInputPort('act_input')
    baitset_input = ObjectInputPort('baitset_input')
    success = ObjectOutputPort('success')

    def begin(self):
        self.baitsets = list()
        self.act_list = list()
        self.ranking = list()

    def process(self, data, port):
        if port is 'act_input':
            self.act_list = data
            if self.args.url != '':
                self.dataset_infos = self.add_dataset()
            
        if port is 'baitset_input':
            self.baitsets.append(data)

        if len(self.act_list) > 0 :
            while len(self.baitsets) > 0 :
                if self.args.url != '':
                    self.success.emit((self.act_list, self.baitsets.pop(), self.ranking, self.dataset_infos))
                else:
                    self.success.emit((self.act_list, self.baitsets.pop(), self.ranking))

    def add_dataset(self):
        url = self.args.url + "/datasets/"
        act_mol_idx = {}
        dataset = None
        parameters = {}

        self.dataset = tempfile.NamedTemporaryFile(suffix='.oeb', mode='wb', delete=False)
        with oechem.oemolostream(self.dataset.name) as ofs:
            for idx, mol in enumerate(self.act_list):
                act_mol_idx[mol.GetTitle()] = idx
                oechem.OEWriteMolecule(ofs, mol)
        self.dataset.flush()

        dataset = open(self.dataset.name, 'rb')
        parameters["dataset"] = (self.dataset.name, dataset, 'application/octet-stream')

        letters = "abcdefghijklmnopqrstuvwxyz123456789"
        uniqueness = ""
        for i in range(0, 8):
            uniqueness += letters[random.randint(0, len(letters)-1)]
        parameters["name"] = 'dataset of active molecules ' + uniqueness

        multipart_data = MultipartEncoder(
            fields=parameters
        )

        response = requests.post(
            url,
            data=multipart_data,
            headers={"content-type": multipart_data.content_type}
        )

        if dataset is not None:
            dataset.close()
        os.remove(self.dataset.name)
 
        data = response.json()
        dataset_infos = (data["id"], act_mol_idx)
        return dataset_infos

class ParallelTreeFPRanking(ParallelComputeCube):
    """
    A compute Cube that receives a Molecule and a list of Fingerprints with a baitset of indices
    and returns the max Similarity value of the Molecule against the Fingerprints
    """

    classification = [["Compute", "Fingerprint", "Similarity"]]

    url = parameter.StringParameter('url', default="http://10.0.62.124:8081", help_text="Url of the FastFingerPrint Server for the request")

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    data_input = ObjectInputPort('data_input')
    success = ObjectOutputPort('success')


    def begin(self):
        pass

    def process(self, data, port):

        self.act_list = data[0]
        self.baitset = data[1]
        self.ranking = data[2]
        for idx in self.baitset[1]:
            smiles = oechem.OEMolToSmiles(self.act_list[idx])
            safe_smiles = parse.quote(smiles)
            url = "%s/%s/hitlist?smiles=%s&oformat=csv&maxhits=%d" %(self.args.url, 'tree_db', safe_smiles, self.args.topn) 
            response = requests.get( url )
            hitlist = response.content.decode().split('\n')
            sys.stdout.flush()
            hitlist.pop(0)
            hitlist.pop()
            cur_rank = list()
            for mol in hitlist:
                cur_mol = mol.split(',')
                cur_rank.append((cur_mol[1], float(cur_mol[2]), self.baitset[0], False))
            if len(self.ranking) == 0:
                self.ranking = cur_rank
            else:
                self.merge_ranking(cur_rank)

        self.success.emit((self.act_list, self.baitset, self.ranking))

    def merge_ranking(self, ranking):
        merged_list = list()
        i = 0
        j = 0
        count = 0
        id_set = set()
        while i < len(self.ranking):
            while j < len(ranking) and ranking[j][1] > self.ranking[i][1]:
                if ranking[j][0] not in id_set: 
                    if count < self.args.topn or ranking[j][1] == merged_list[count-1][1]:
                        merged_list.append(ranking[j])
                        count += 1
                        id_set.add(ranking[j][0])
                        j += 1
                    else:
                        break
                else:
                    j += 1

            if self.ranking[i][0] not in id_set: 
                if self.ranking[i] not in id_set and (count < self.args.topn or self.ranking[i][1] == merged_list[count-1][1]):
                    merged_list.append(self.ranking[i])  
                    count += 1
                    id_set.add(self.ranking[i][0])
                    i += 1
                else:
                    break
            else:
                i += 1

        while j < len(ranking):
            if ranking[j][0] not in id_set: 
                if ranking[j] not in id_set and (count < self.args.topn or ranking[j][1] == merged_list[count-1][1]):
                    merged_list.append(ranking[j])
                    count += 1
                    id_set.add(ranking[j][0])
                    j += 1
                else:
                    break
            else:
                j += 1

        self.ranking = merged_list

    def update_ranking(self, mol, max_tanimoto, ka_tag):
        index = 0
        if len(self.ranking) >= self.args.topn and max_tanimoto < self.ranking[len(self.ranking)-1][1]:
            pass
        else:    
            for top_mol in self.ranking:
                if max_tanimoto < top_mol[1]:
                    index = self.ranking.index(top_mol) + 1
                else:
                    break

            upper = self.ranking[:index]
            lower = self.ranking[index:]
            self.ranking = upper + [(mol.GetTitle(), max_tanimoto, self.baitset[0], ka_tag)] + lower

            i = self.args.topn - 1
            while i < len(self.ranking) - 1:
                if self.ranking[i][1] != self.ranking[i + 1][1]:
                    self.ranking = self.ranking[:i + 1]

                    break
                else:
                    i += 1

class ParallelTreeFPInsertKA(ParallelComputeCube):
    """
    """

    classification = [["ParallelCompute"]]

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    numbits = parameter.IntegerParameter('numbits', default=4096,
                                    help_text="size of the fingerprints")

    data_input = ObjectInputPort('data_input')
    success = ObjectOutputPort('success')

    def process(self, data, port):
        
        self.act_list = data[0]
        self.baitset = data[1]
        self.ranking = data[2]
        self.fptype = oegraphsim.OEFPType_Tree 
        self.fp_list = list()

        self.calculate_fp()
        self.insert_known_actives()

        self.success.emit((self.act_list, self.baitset, self.ranking)) 

    def calculate_fp(self):
        minbond = 0
        maxbond = 4
        atype = AtmNum|Arom|Chiral|FCharge|HvyDeg|Hyb|EqHalo
        btype = Order|ChiralB

        for mol in self.act_list:
            fp = oegraphsim.OEFingerPrint()
            oegraphsim.OEMakeTreeFP(fp, mol, self.args.numbits, minbond, maxbond, atype, btype)
            #oegraphsim.OEMakeFP(fp, mol, self.fptype)

            self.fp_list.append(fp)

    def insert_known_actives(self):

        c = 0
        for idx in self.baitset[1]:
            while c < idx:
                ka_fp = self.fp_list[c]
                simval = self.calc_sim_val(ka_fp)
                self.update_ranking(self.act_list[c], simval, True)

                c += 1
            c += 1
        while c < len(self.act_list):
            ka_fp = self.fp_list[c]
            simval = self.calc_sim_val(ka_fp)
            self.update_ranking(self.act_list[c], simval, True)
            c += 1

    def calc_sim_val(self, fp):
        maxval = 0
        for idx in self.baitset[1]:
            tanimoto = oechem.OETanimoto(fp, self.fp_list[idx])
            if tanimoto > maxval:
                maxval = tanimoto
        return maxval

    def update_ranking(self, mol, max_tanimoto, ka_tag):
        index = 0
        if len(self.ranking) >= self.args.topn and max_tanimoto < self.ranking[len(self.ranking)-1][1]:
            pass
        else:    
            for top_mol in self.ranking:
                if max_tanimoto < top_mol[1]:
                    index = self.ranking.index(top_mol) + 1
                else:
                    break

            upper = self.ranking[:index]
            lower = self.ranking[index:]
            self.ranking = upper + [(mol.GetTitle(), max_tanimoto, self.baitset[0], ka_tag)] + lower

            i = self.args.topn - 1
            while i < len(self.ranking) - 1:
                if self.ranking[i][1] != self.ranking[i + 1][1]:
                    self.ranking = self.ranking[:i + 1]

                    break
                else:
                    i += 1

    def end(self):
        pass

class ParallelPathFPRanking(ParallelComputeCube):
    """
    A compute Cube that receives a Molecule and a list of Fingerprints with a baitset of indices
    and returns the max Similarity value of the Molecule against the Fingerprints
    """

    classification = [["Compute", "Fingerprint", "Similarity"]]

    url = parameter.StringParameter('url', default="http://10.0.62.124:8081", help_text="Url of the FastFingerPrint Server for the request")

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    data_input = ObjectInputPort('data_input')
    success = ObjectOutputPort('success')


    def begin(self):
        pass

    def process(self, data, port):

        self.act_list = data[0]
        self.baitset = data[1]
        self.ranking = data[2]
        for idx in self.baitset[1]:
            smiles = oechem.OEMolToSmiles(self.act_list[idx])
            safe_smiles = parse.quote(smiles)
            url = "%s/%s/hitlist?smiles=%s&oformat=csv&maxhits=%d" %(self.args.url, 'path_db', safe_smiles, self.args.topn) 
            response = requests.get( url )
            hitlist = response.content.decode().split('\n')
            hitlist.pop(0)
            hitlist.pop()
            cur_rank = list()
            for mol in hitlist:
                cur_mol = mol.split(',')
                cur_rank.append((cur_mol[1], float(cur_mol[2]), self.baitset[0], False))
            if len(self.ranking) == 0:
                self.ranking = cur_rank
            else:
                self.merge_ranking(cur_rank)

        self.success.emit((self.act_list, self.baitset, self.ranking))

    def merge_ranking(self, ranking):
        merged_list = list()
        i = 0
        j = 0
        count = 0
        id_set = set()
        while i < len(self.ranking):
            while j < len(ranking) and ranking[j][1] > self.ranking[i][1]:
                if ranking[j][0] not in id_set: 
                    if count < self.args.topn or ranking[j][1] == merged_list[count-1][1]:
                        merged_list.append(ranking[j])
                        count += 1
                        id_set.add(ranking[j][0])
                        j += 1
                    else:
                        break
                else:
                    j += 1

            if self.ranking[i][0] not in id_set: 
                if self.ranking[i] not in id_set and (count < self.args.topn or self.ranking[i][1] == merged_list[count-1][1]):
                    merged_list.append(self.ranking[i])  
                    count += 1
                    id_set.add(self.ranking[i][0])
                    i += 1
                else:
                    break
            else:
                i += 1

        while j < len(ranking):
            if ranking[j][0] not in id_set: 
                if ranking[j] not in id_set and (count < self.args.topn or ranking[j][1] == merged_list[count-1][1]):
                    merged_list.append(ranking[j])
                    count += 1
                    id_set.add(ranking[j][0])
                    j += 1
                else:
                    break
            else:
                j += 1

        self.ranking = merged_list

    def update_ranking(self, mol, max_tanimoto, ka_tag):
        index = 0
        if len(self.ranking) >= self.args.topn and max_tanimoto < self.ranking[len(self.ranking)-1][1]:
            pass
        else:    
            for top_mol in self.ranking:
                if max_tanimoto < top_mol[1]:
                    index = self.ranking.index(top_mol) + 1
                else:
                    break

            upper = self.ranking[:index]
            lower = self.ranking[index:]
            self.ranking = upper + [(mol.GetTitle(), max_tanimoto, self.baitset[0], ka_tag)] + lower

            i = self.args.topn - 1
            while i < len(self.ranking) - 1:
                if self.ranking[i][1] != self.ranking[i + 1][1]:
                    self.ranking = self.ranking[:i + 1]

                    break
                else:
                    i += 1

class ParallelPathFPInsertKA(ParallelComputeCube):
    """
    """

    classification = [["ParallelCompute"]]

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    numbits = parameter.IntegerParameter('numbits', default=4096,
                                    help_text="size of the fingerprints")

    data_input = ObjectInputPort('data_input')
    success = ObjectOutputPort('success')

    def process(self, data, port):
        
        self.act_list = data[0]
        self.baitset = data[1]
        self.ranking = data[2]
        self.fptype = oegraphsim.OEFPType_Path
        self.fp_list = list()

        self.calculate_fp()
        self.insert_known_actives()

        self.success.emit((self.act_list, self.baitset, self.ranking)) 

    def calculate_fp(self):
        minbond = 0
        maxbond = 4
        atype = AtmNum|Arom|Chiral|FCharge|HvyDeg|Hyb|EqHalo
        btype = Order|ChiralB

        for mol in self.act_list:
            fp = oegraphsim.OEFingerPrint()
            oegraphsim.OEMakePathFP(fp, mol, self.args.numbits, minbond, maxbond, atype, btype)
            #oegraphsim.OEMakeFP(fp, mol, self.fptype)

            self.fp_list.append(fp)

    def insert_known_actives(self):

        c = 0
        for idx in self.baitset[1]:
            while c < idx:
                ka_fp = self.fp_list[c]
                simval = self.calc_sim_val(ka_fp)
                self.update_ranking(self.act_list[c], simval, True)

                c += 1
            c += 1
        while c < len(self.act_list):
            ka_fp = self.fp_list[c]
            simval = self.calc_sim_val(ka_fp)
            self.update_ranking(self.act_list[c], simval, True)
            c += 1

    def calc_sim_val(self, fp):
        maxval = 0
        for idx in self.baitset[1]:
            tanimoto = oechem.OETanimoto(fp, self.fp_list[idx])
            if tanimoto > maxval:
                maxval = tanimoto
        return maxval

    def update_ranking(self, mol, max_tanimoto, ka_tag):
        index = 0
        if len(self.ranking) >= self.args.topn and max_tanimoto < self.ranking[len(self.ranking)-1][1]:
            pass
        else:    
            for top_mol in self.ranking:
                if max_tanimoto < top_mol[1]:
                    index = self.ranking.index(top_mol) + 1
                else:
                    break

            upper = self.ranking[:index]
            lower = self.ranking[index:]
            self.ranking = upper + [(mol.GetTitle(), max_tanimoto, self.baitset[0], ka_tag)] + lower

            i = self.args.topn - 1
            while i < len(self.ranking) - 1:
                if self.ranking[i][1] != self.ranking[i + 1][1]:
                    self.ranking = self.ranking[:i + 1]

                    break
                else:
                    i += 1

    def end(self):
        pass

class ParallelCircularFPRanking(ParallelComputeCube):
    """
    A compute Cube that receives a Molecule and a list of Fingerprints with a baitset of indices
    and returns the max Similarity value of the Molecule against the Fingerprints
    """

    classification = [["Compute", "Fingerprint", "Similarity"]]

    url = parameter.StringParameter('url', default="http://10.0.62.124:8081", help_text="Url of the FastFingerPrint Server for the request")

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    data_input = ObjectInputPort('data_input')
    success = ObjectOutputPort('success')


    def begin(self):
        pass

    def process(self, data, port):

        self.act_list = data[0]
        self.baitset = data[1]
        self.ranking = data[2]
        for idx in self.baitset[1]:
            smiles = oechem.OEMolToSmiles(self.act_list[idx])
            safe_smiles = parse.quote(smiles)
            url = "%s/%s/hitlist?smiles=%s&oformat=csv&maxhits=%d" %(self.args.url, 'circular_db', safe_smiles, self.args.topn) 
            response = requests.get( url )
            hitlist = response.content.decode().split('\n')
            sys.stdout.flush()
            hitlist.pop(0)
            hitlist.pop()
            cur_rank = list()
            for mol in hitlist:
                cur_mol = mol.split(',')
                cur_rank.append((cur_mol[1], float(cur_mol[2]), self.baitset[0], False))
            if len(self.ranking) == 0:
                self.ranking = cur_rank
            else:
                self.merge_ranking(cur_rank)

        self.success.emit((self.act_list, self.baitset, self.ranking))

    def merge_ranking(self, ranking):
        merged_list = list()
        i = 0
        j = 0
        count = 0
        id_set = set()
        while i < len(self.ranking):
            while j < len(ranking) and ranking[j][1] > self.ranking[i][1]:
                if ranking[j][0] not in id_set: 
                    if count < self.args.topn or ranking[j][1] == merged_list[count-1][1]:
                        merged_list.append(ranking[j])
                        count += 1
                        id_set.add(ranking[j][0])
                        j += 1
                    else:
                        break
                else:
                    j += 1

            if self.ranking[i][0] not in id_set: 
                if self.ranking[i] not in id_set and (count < self.args.topn or self.ranking[i][1] == merged_list[count-1][1]):
                    merged_list.append(self.ranking[i])  
                    count += 1
                    id_set.add(self.ranking[i][0])
                    i += 1
                else:
                    break
            else:
                i += 1

        while j < len(ranking):
            if ranking[j][0] not in id_set: 
                if ranking[j] not in id_set and (count < self.args.topn or ranking[j][1] == merged_list[count-1][1]):
                    merged_list.append(ranking[j])
                    count += 1
                    id_set.add(ranking[j][0])
                    j += 1
                else:
                    break
            else:
                j += 1

        self.ranking = merged_list

    def update_ranking(self, mol, max_tanimoto, ka_tag):
        index = 0
        if len(self.ranking) >= self.args.topn and max_tanimoto < self.ranking[len(self.ranking)-1][1]:
            pass
        else:    
            for top_mol in self.ranking:
                if max_tanimoto < top_mol[1]:
                    index = self.ranking.index(top_mol) + 1
                else:
                    break

            upper = self.ranking[:index]
            lower = self.ranking[index:]
            self.ranking = upper + [(mol.GetTitle(), max_tanimoto, self.baitset[0], ka_tag)] + lower

            i = self.args.topn - 1
            while i < len(self.ranking) - 1:
                if self.ranking[i][1] != self.ranking[i + 1][1]:
                    self.ranking = self.ranking[:i + 1]

                    break
                else:
                    i += 1

class ParallelCircularFPInsertKA(ParallelComputeCube):
    """
    """

    classification = [["ParallelCompute"]]

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    numbits = parameter.IntegerParameter('numbits', default=4096,
                                    help_text="size of the fingerprints")

    data_input = ObjectInputPort('data_input')
    success = ObjectOutputPort('success')

    def process(self, data, port):
        
        self.act_list = data[0]
        self.baitset = data[1]
        self.ranking = data[2]
        self.fptype = oegraphsim.OEFPType_Circular
        self.fp_list = list()

        self.calculate_fp()
        self.insert_known_actives()

        self.success.emit((self.act_list, self.baitset, self.ranking)) 

    def calculate_fp(self):
        minbond = 0
        maxbond = 4
        atype = AtmNum|Arom|Chiral|FCharge|HvyDeg|Hyb|EqHalo
        btype = Order|ChiralB
        
        for mol in self.act_list:
            fp = oegraphsim.OEFingerPrint()
            oegraphsim.OEMakeCircularFP(fp, mol, self.args.numbits, minbond, maxbond, atype, btype)
            #oegraphsim.OEMakeFP(fp, mol, self.fptype)

            self.fp_list.append(fp)

    def insert_known_actives(self):

        c = 0
        for idx in self.baitset[1]:
            while c < idx:
                ka_fp = self.fp_list[c]
                simval = self.calc_sim_val(ka_fp)
                self.update_ranking(self.act_list[c], simval, True)

                c += 1
            c += 1
        while c < len(self.act_list):
            ka_fp = self.fp_list[c]
            simval = self.calc_sim_val(ka_fp)
            self.update_ranking(self.act_list[c], simval, True)
            c += 1

    def calc_sim_val(self, fp):
        maxval = 0
        for idx in self.baitset[1]:
            tanimoto = oechem.OETanimoto(fp, self.fp_list[idx])
            if tanimoto > maxval:
                maxval = tanimoto
        return maxval

    def update_ranking(self, mol, max_tanimoto, ka_tag):
        index = 0
        if len(self.ranking) >= self.args.topn and max_tanimoto < self.ranking[len(self.ranking)-1][1]:
            pass
        else:    
            for top_mol in self.ranking:
                if max_tanimoto < top_mol[1]:
                    index = self.ranking.index(top_mol) + 1
                else:
                    break

            upper = self.ranking[:index]
            lower = self.ranking[index:]
            self.ranking = upper + [(mol.GetTitle(), max_tanimoto, self.baitset[0], ka_tag)] + lower

            i = self.args.topn - 1
            while i < len(self.ranking) - 1:
                if self.ranking[i][1] != self.ranking[i + 1][1]:
                    self.ranking = self.ranking[:i + 1]

                    break
                else:
                    i += 1

    def end(self):
        pass

class ParallelFastROCSRanking(ParallelComputeCube):
    """
    A compute Cube that receives a Molecule a baitset of indices and a FastROCSServer address
    and returns the ranking of the Server Molecules against the query
    """

    classification = [["Compute", "FastROCS", "Similarity"]]

    url = parameter.StringParameter('url', default="http://10.0.61.25:4711", help_text="Url of the FastROCS Server for the request")

    dataset_name = parameter.StringParameter('dataset_name', default="screening_database", help_text="Name of the screening database")

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    wait = parameter.IntegerParameter('wait', default=1,
                                    help_text="waiting time before trying a new request")

    data_input = ObjectInputPort('data_input')
    success = ObjectOutputPort('success')

    def begin(self):
        pass

    def process(self, data, port):

        self.act_list = data[0]
        self.baitset = data[1]
        self.ranking = data[2]
        self.dataset_infos = data[3]

        self.log.info("start ranking baitset number {}".format(self.baitset[0]))

        url = self.args.url + "/datasets/"
        response = requests.get(url)
        data = response.json()
        datasets = data["datasets"]
        for dataset in datasets:
            if dataset["name"] == self.args.dataset_name:
                self.dataset_identifier = int(dataset["id"])

        count = 0
        results = self.run_query()

        cur_ranking_dict = self.create_rankings(results)

        for cur_rank in cur_ranking_dict.values():
            if len(self.ranking) == 0:
                self.ranking = cur_rank
            else:
                self.merge_ranking(cur_rank)
            count += 1
            self.log.info("Baitset " + str(self.baitset[0]) + " : " + str(count) +" requests processed")

        self.log.info("Emitting ranking baitset " + str(self.baitset[0]))
        self.success.emit((self.act_list, self.baitset, self.ranking, self.dataset_infos)) 

    def run_query(self):
        url = self.args.url + "/queries/"

        self.query = tempfile.NamedTemporaryFile(suffix='.oeb', mode='wb', delete=False)  
        with oechem.oemolostream(self.query.name) as ofs:
            for idx in self.baitset[1]:
                oechem.OEWriteMolecule(ofs, self.act_list[idx])
        self.query.flush()

        parameters = {}
        parameters["num_hits"] = self.args.topn
        parameters["dataset_identifier"] = self.dataset_identifier
        parameters["cutoff"] = 0.0
        parameters["ScaffoldCutoff"] = 1.0
        parameters["SimFunc"] = 'Tanimoto'
        parameters["SimType"] = 'Combo'
        parameters["shape_only"] = False

        with open(self.query.name, "rb") as query_file:
            response = requests.post(
                url,
                files={"query": query_file},
                data=parameters
            )
            if not response.ok:
                print(response.json()["error"])
                return

        os.remove(self.query.name)
        data = response.json()
        query_id = data["id"]

        status_url = url + "{}/".format(query_id)
        tries = 0
        while True:
            time.sleep(self.args.wait*tries)
            response = requests.get(status_url)
            if not response.ok:
                print(response.json()["error"])
                return
            results = response.json()
            status = results["status"]
            total = int(status["total"])
            current = int(status["current"])

            if status["job"] == "FAILED":
                print("Query {0:d} failed".format(query_id))
                return

            if(total != 0 and current >= total and
               status["job"] == "COMPLETED"):
                break

        response = requests.get(self.args.url + results["results"])

        requests.delete(status_url)
        return response

    def create_rankings(self, results_data):
        cur_ranking_dict = dict()
        with tempfile.NamedTemporaryFile(suffix='.oeb', mode='wb', delete=False) as temp:
            temp.write(results_data.content)
            temp.flush()
            with oechem.oemolistream(temp.name) as results:
                for mol in results.GetOEGraphMols():
                    if oechem.OEGetSDData(mol, 'QueryMol') in cur_ranking_dict.keys():
                        if mol.GetTitle() not in [mol[0] for mol in cur_ranking_dict[oechem.OEGetSDData(mol, 'QueryMol')]]:
                            cur_ranking_dict[oechem.OEGetSDData(mol, 'QueryMol')].append((mol.GetTitle(), float(oechem.OEGetSDData(mol, 'TanimotoCombo')), self.baitset[0], False))
                    else:
                        cur_rank = list()
                        cur_rank.append((mol.GetTitle(), float(oechem.OEGetSDData(mol, 'TanimotoCombo')), self.baitset[0], False))
                        cur_ranking_dict[oechem.OEGetSDData(mol, 'QueryMol')] = cur_rank

            os.remove(temp.name)
            return cur_ranking_dict

    def merge_ranking(self, ranking):
        merged_list = list()
        i = 0
        j = 0
        count = 0
        id_set = set()
        while i < len(self.ranking):
            while j < len(ranking) and ranking[j][1] > self.ranking[i][1]:
                if ranking[j][1] not in id_set: 
                    if count < self.args.topn or ranking[j][1] == merged_list[count-1][1]:
                        merged_list.append(ranking[j])
                        count += 1
                        id_set.add(ranking[j][0])
                        j += 1
                    else:
                        break
                else:
                    j += 1

            if self.ranking[i][0] not in id_set: 
                if self.ranking[i] not in id_set and (count < self.args.topn or self.ranking[i][1] == merged_list[count-1][1]):
                    merged_list.append(self.ranking[i])  
                    count += 1
                    id_set.add(self.ranking[i][0])
                    i += 1
                else:
                    break
            else:
                i += 1

        while j < len(ranking):
            if ranking[j][0] not in id_set: 
                if ranking[j] not in id_set and (count < self.args.topn or ranking[j][1] == merged_list[count-1][1]):
                    merged_list.append(ranking[j])
                    count += 1
                    id_set.add(ranking[j][0])
                    j += 1
                else:
                    break
            else:
                j += 1

        self.ranking = merged_list

class ParallelInsertKARestfulROCS(ParallelComputeCube):
    """
    """

    classification = [["ParallelCompute"]]

    url = parameter.StringParameter('url', default="http://10.0.61.25:4711", help_text="Url of the Restful FastROCS Server for the request")

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    data_input = ObjectInputPort('data_input')
    success = ObjectOutputPort('success')

    def process(self, data, port):
        
        self.act_list = data[0]
        self.baitset = data[1]
        self.ranking = data[2]
        self.dataset_infos = data[3]
        self.log.info("processing KA for baitset : " + str(self.baitset[0]))

        self.dataset_identifier = self.dataset_infos[0]

        results = self.run_query()

        self.create_cur_scores(results)
        for tanimoto, mol in self.cur_scores.values():
                self.update_ranking(mol, tanimoto, True)

        self.success.emit((self.act_list, self.baitset, self.ranking, self.dataset_infos[0]))

    def run_query(self):
        url = self.args.url + "/queries/"

        self.query = tempfile.NamedTemporaryFile(suffix='.oeb', mode='wb', delete=False)  
        with oechem.oemolostream(self.query.name) as ofs:
            for idx in self.baitset[1]:
                oechem.OEWriteMolecule(ofs, self.act_list[idx])
        self.query.flush()

        parameters = {}
        parameters["num_hits"] = self.args.topn
        parameters["dataset_identifier"] = self.dataset_identifier
        parameters["cutoff"] = 0.0
        parameters["ScaffoldCutoff"] = 1.0
        parameters["SimFunc"] = 'Tanimoto'
        parameters["SimType"] = 'Combo'
        parameters["shape_only"] = False
        parameters["oformat"] = 'oeb'

        with open(self.query.name, "rb") as query_file:
            response = requests.post(
                url,
                files={"query": query_file},
                data=parameters
            )
            if not response.ok:
                print(response.json()["error"])
                return

        os.remove(self.query.name)
        data = response.json()
        query_id = data["id"]

        status_url = url + "{}/".format(query_id)
        tries = 0
        while True:
            time.sleep(tries)
            tries += 1
            response = requests.get(status_url)
            if not response.ok:
                print(response.json()["error"])
                return
            results = response.json()
            status = results["status"]
            total = int(status["total"])
            current = int(status["current"])

            if status["job"] == "FAILED":
                print("Query {0:d} failed".format(query_id))
                return

            if(total != 0 and current >= total and
               status["job"] == "COMPLETED"):
                break

        response = requests.get(self.args.url + results["results"])

        requests.delete(status_url)
        return response

    def create_cur_scores(self, results_data):
        self.cur_scores = {}
#        bait=False
#        baitname=''
        with tempfile.NamedTemporaryFile(suffix='.oeb', mode='wb', delete=False) as temp:
            temp.write(results_data.content)
            temp.flush()
            with oechem.oemolistream(temp.name) as results:
                for mol in results.GetOEGraphMols():
                    if self.dataset_infos[1][mol.GetTitle()] not in self.baitset[1]:
                        tanimoto_combo = float(oechem.OEGetSDData(mol, "TanimotoCombo"))
                        if mol.GetTitle() in self.cur_scores.keys():
                            if self.cur_scores[mol.GetTitle()][0] < tanimoto_combo:
                                self.cur_scores[mol.GetTitle()] = (tanimoto_combo, mol.CreateCopy())
                        else:
                            self.cur_scores[mol.GetTitle()] = (tanimoto_combo, mol.CreateCopy())
                    #Test/#
#                    else:
#                        if bait!=True or mol.GetTitle() == baitname:
#                            baitname = mol.GetTitle()
#                            print(self.baitset[0], mol.GetTitle())
#                            tanimoto_combo = float(oechem.OEGetSDData(mol, "TanimotoCombo"))
#                            if mol.GetTitle() in self.cur_scores.keys():
#                                if self.cur_scores[mol.GetTitle()][0] < tanimoto_combo:
#                                    self.cur_scores[mol.GetTitle()] = (tanimoto_combo, mol.CreateCopy())
#                            else:
#                                self.cur_scores[mol.GetTitle()] = (tanimoto_combo, mol.CreateCopy())
#                            bait=True
                     #/Test#
                            
            os.remove(temp.name)

    def update_ranking(self, mol, max_tanimoto, ka_tag):
        index = 0
        if len(self.ranking) >= self.args.topn and max_tanimoto < self.ranking[len(self.ranking)-1][1]:
            pass
        else:    
            for top_mol in self.ranking:
                if max_tanimoto < top_mol[1]:
                    index = self.ranking.index(top_mol) + 1
                else:
                    break

            upper = self.ranking[:index]
            lower = self.ranking[index:]
            self.ranking = upper + [(mol.GetTitle(), max_tanimoto, self.baitset[0], ka_tag)] + lower

            i = self.args.topn - 1
            while i < len(self.ranking) - 1:
                if self.ranking[i][1] != self.ranking[i + 1][1]:
                    self.ranking = self.ranking[:i + 1]

                    break
                else:
                    i += 1

    def end(self):
        pass

class AccumulateRankings(ComputeCube):
    """
    A compute Cube that receives rankings and assemble them in a list
    """

    classification = [["Compute", "Accumulator"]]

    url = parameter.StringParameter('url', default="http://10.0.1.22:4242", help_text="Url of the Restful FastROCS Server for the request")

    tree_fpintake = ObjectInputPort('tree_fpintake')
    path_fpintake = ObjectInputPort('path_fpintake')
    circular_fpintake = ObjectInputPort('circular_fpintake')
    rocsintake = ObjectInputPort('rocsintake')
    success = ObjectOutputPort('success')

    def begin(self):
        self.tree_fp_ranking_list = list()
        self.path_fp_ranking_list = list()
        self.circular_fp_ranking_list = list()
        self.rocs_ranking_list = list()

    def process(self, data, port):
        if port == 'tree_fpintake':
            self.tree_fp_ranking_list.append(data[2])
            self.nb_ka = len(data[0])-len(data[1][1])
        if port == 'path_fpintake':
            self.path_fp_ranking_list.append(data[2])
            self.nb_ka = len(data[0])-len(data[1][1])
        if port == 'circular_fpintake':
            self.circular_fp_ranking_list.append(data[2])
            self.nb_ka = len(data[0])-len(data[1][1])
        if port == 'rocsintake':
            self.rocs_ranking_list.append(data[2])
            self.nb_ka = len(data[0])-len(data[1][1])
            self.dataset_id = data[3]

    def end(self):
        url = self.args.url + '/datasets/'
        if self.dataset_id != None:
            response = requests.delete(url + str(self.dataset_id) + '/')
        self.success.emit((self.tree_fp_ranking_list, self.nb_ka, 'Tree_FP')) 
        self.success.emit((self.path_fp_ranking_list, self.nb_ka, 'Path_FP')) 
        self.success.emit((self.circular_fp_ranking_list, self.nb_ka, 'Circular_FP')) 
        self.success.emit((self.rocs_ranking_list, self.nb_ka, 'FastROCS'))

class AccumulateRankingsFP(ComputeCube):
    """
    A compute Cube that receives rankings and assemble them in a list
    """

    classification = [["Compute", "Accumulator"]]

    tree_fpintake = ObjectInputPort('tree_fpintake')
    path_fpintake = ObjectInputPort('path_fpintake')
    circular_fpintake = ObjectInputPort('circular_fpintake')
    success = ObjectOutputPort('success')

    def begin(self):
        self.tree_fp_ranking_list = list()
        self.path_fp_ranking_list = list()
        self.circular_fp_ranking_list = list()

    def process(self, data, port):
        if port == 'tree_fpintake':
            self.tree_fp_ranking_list.append(data[2])
            self.nb_ka = len(data[0])-len(data[1][1])
        if port == 'path_fpintake':
            self.path_fp_ranking_list.append(data[2])
            self.nb_ka = len(data[0])-len(data[1][1])
        if port == 'circular_fpintake':
            self.circular_fp_ranking_list.append(data[2])
            self.nb_ka = len(data[0])-len(data[1][1])

    def end(self):
        self.success.emit((self.tree_fp_ranking_list, self.nb_ka, 'Tree_FP')) 
        self.success.emit((self.path_fp_ranking_list, self.nb_ka, 'Path_FP')) 
        self.success.emit((self.circular_fp_ranking_list, self.nb_ka, 'Circular_FP')) 

class AccumulateRankingsTree(ComputeCube):
    """
    A compute Cube that receives rankings and assemble them in a list
    """

    classification = [["Compute", "Accumulator"]]

    tree_fpintake = ObjectInputPort('tree_fpintake')
    success = ObjectOutputPort('success')

    def begin(self):
        self.tree_fp_ranking_list = list()

    def process(self, data, port):
        if port == 'tree_fpintake':
            self.tree_fp_ranking_list.append(data[2])
            self.nb_ka = len(data[0])-len(data[1][1])

    def end(self):
        self.success.emit((self.tree_fp_ranking_list, self.nb_ka, 'Tree_FP')) 

class AccumulateRankingsPath(ComputeCube):
    """
    A compute Cube that receives rankings and assemble them in a list
    """

    classification = [["Compute", "Accumulator"]]

    path_fpintake = ObjectInputPort('path_fpintake')
    success = ObjectOutputPort('success')

    def begin(self):
        self.path_fp_ranking_list = list()

    def process(self, data, port):
        if port == 'path_fpintake':
            self.path_fp_ranking_list.append(data[2])
            self.nb_ka = len(data[0])-len(data[1][1])

    def end(self):
        self.success.emit((self.path_fp_ranking_list, self.nb_ka, 'Path_FP')) 

class AccumulateRankingsCircular(ComputeCube):
    """
    A compute Cube that receives rankings and assemble them in a list
    """

    classification = [["Compute", "Accumulator"]]

    circular_fpintake = ObjectInputPort('circular_fpintake')
    success = ObjectOutputPort('success')

    def begin(self):
        self.circular_fp_ranking_list = list()

    def process(self, data, port):
        if port == 'circular_fpintake':
            self.circular_fp_ranking_list.append(data[2])
            self.nb_ka = len(data[0])-len(data[1][1])

    def end(self):
        self.success.emit((self.circular_fp_ranking_list, self.nb_ka, 'Circular_FP')) 

class AccumulateRankingsROCS(ComputeCube):
    """
    A compute Cube that receives rankings and assemble them in a list
    """

    classification = [["Compute", "Accumulator"]]

    url = parameter.StringParameter('url', default="http://10.0.1.22:4242", help_text="Url of the Restful FastROCS Server for the request")

    rocsintake = ObjectInputPort('rocsintake')
    success = ObjectOutputPort('success')

    def begin(self):
        self.rocs_ranking_list = list()

    def process(self, data, port):
        if port == 'rocsintake':
            self.rocs_ranking_list.append(data[2])
            self.nb_ka = len(data[0])-len(data[1][1])
            self.dataset_id = data[3]

    def end(self):
        url = self.args.url + '/datasets/'
        if self.dataset_id != None:
            response = requests.delete(url + str(self.dataset_id) + '/')
        self.success.emit((self.rocs_ranking_list, self.nb_ka, 'FastROCS'))

class AnalyseRankings(ComputeCube):
    """

    """

    classification = [["Compute", "Analysis"]]

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    intake = ObjectInputPort('intake')
    success = ObjectOutputPort('success')

    def process(self, data, port):
        self.ranking_list = data[0]
        self.nb_ka = data[1]
        self.method = data[2]
        self.datas = self.ranking_analysis()

        self.success.emit((self.datas, self.method))

    def ranking_analysis(self):
        results = pd.DataFrame()
        for ranking in self.ranking_list:
            set_results = pd.DataFrame(columns = ['RR', 'HR'])
            count = 0
            count_ka = 0
            for row, mol in enumerate(ranking):
                count += 1
                if mol[3] == 1:
                    count_ka += 1
                rr = 100 * count_ka/self.nb_ka
                hr = 100 * count_ka/count
                set_results.loc[row] = [rr, hr]
            results = pd.concat([results, set_results])
        
        datas = pd.DataFrame()

        if self.method == 'Tree_FP':
            name = 'FP_tree'
        if self.method == 'Path_FP':
            name = 'FP_path'
        if self.method == 'Circular_FP':
            name = 'FP_circular'
        elif self.method == 'FastROCS':
            name = 'FR'

        datas['Average RR ' + name] = results.groupby(results.index)['RR'].mean()
        datas['Average HR ' + name] = results.groupby(results.index)['HR'].mean()

        #datas['Minimum RR ' + name] = results.groupby(results.index)['RR'].min()
        #datas['Minimum HR ' + name] = results.groupby(results.index)['HR'].min()

        #datas['Maximum RR ' + name] = results.groupby(results.index)['RR'].max()
        #datas['Maximum HR ' + name] = results.groupby(results.index)['HR'].max()
        datas = datas.head(self.args.topn)

        return datas

