
import json,requests
import urllib.parse as parse

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from openeye import oechem
from openeye import oegraphsim
from openeye import oeomega
from openeye import oemolprop

import time
import random

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
        self.log.info("Generator Cube runnig")

    def process(self, data, port):
        pass

    def end(self):
        total = 100
        nb_index = 33 
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


class CalculateFPCube(ComputeCube):
    """
    A compute Cube that reads a list of Molecules and returns a list of Fingerprints
    """

    classification = [["Compute", "Fingerprint"]]

    fptype = parameter.IntegerParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")

    intake = ObjectInputPort('intake')
    success = ObjectOutputPort('success')


    def begin(self):
        self.fp_list = list()

    def process(self, data, port):
        for mol in data:
            fp = oegraphsim.OEFingerPrint()
            oegraphsim.OEMakeFP(fp, mol, self.args.fptype)
#            bitstring = ''
#            for b in range(0, fp.GetSize()):
#                if fp.IsBitOn(b):
#                    bitstring += '1'
#                else:
#                    bitstring += '0'

            self.fp_list.append(fp)

    def end(self):
  #      for fp in self.fp_list:
        self.success.emit(self.fp_list)

class PrepareRanking(ComputeCube):

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

        if port is 'baitset_input':
            self.baitsets.append(data)

        if len(self.act_list) > 0 :
            while len(self.baitsets) > 0 :
                self.success.emit((self.act_list, self.baitsets.pop(), self.ranking))

class ParallelRanking(ParallelComputeCube):
    """
    A compute Cube that receives a Molecule and a list of Fingerprints with a baitset of indices
    and returns the max Similarity value of the Molecule against the Fingerprints
    """

    classification = [["Compute", "Fingerprint", "Similarity"]]

    fptype = parameter.IntegerParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    data_input = ObjectInputPort('data_input')
    success = ObjectOutputPort('success')


    def begin(self):
#        self.max_tanimoto = 0
#        self.fp = None
#        self.fp_list = None
#        self.baitset = None
        pass

    def process(self, data, port):

        self.act_list = data[0]
        self.baitset = data[1]
        self.ranking = data[2]
        baseurl = "http://130.180.63.34:8069"
        fptypes = {102 : 'path', 104 : 'circular', 105 : 'tree'}
        database = fptypes[self.args.fptype] + "_db"
        for idx in self.baitset[1]:
            smiles = oechem.OEMolToSmiles(self.act_list[idx])
            safe_smiles = parse.quote(smiles)
            url = "%s/%s/hitlist?smiles=%s&oformat=csv&maxhits=%d" %(baseurl, database, safe_smiles, self.args.topn) 
            response = requests.get( url )
            hitlist = response.content.decode().split('\n')
            hitlist.pop(0)
            hitlist.pop()
            cur_rank = list()
            for mol in hitlist:
                cur_mol = mol.split(',')
                cur_rank.append((cur_mol[0], cur_mol[1], float(cur_mol[5]), self.baitset[0], False))
            if len(self.ranking) == 0:
                self.ranking = cur_rank
            else:
                self.merge_ranking(cur_rank)

        #if self.fp_list is not None and self.baitset is not None:
        #with oechem.oemolistream(str(self.args.data_in)) as ifs:
        #    for mol in ifs.GetOEMols():
        #        max_tanimoto = 0
        #        fp = oegraphsim.OEFingerPrint()
        #        oegraphsim.OEMakeFP(fp, mol, self.args.fptype)
        #        for idx in self.baitset[1]:
        #            act_fp = self.fp_list[idx]
        #            tanimoto = oegraphsim.OETanimoto(fp, self.fp_list[idx])
        #            if tanimoto > max_tanimoto:
        #                max_tanimoto = tanimoto
        #        self.update_ranking(mol, max_tanimoto, False)

        self.success.emit((self.act_list, self.baitset, self.ranking))

    def merge_ranking(self, ranking):
        merged_list = list()
        i = 0
        j = 0
        count = 0
        id_set = set()
        while i < len(self.ranking):
            while j < len(ranking) and ranking[j][2] > self.ranking[i][2]:
                if ranking[j][1] not in id_set: 
                    if count < self.args.topn or ranking[j][2] == merged_list[count-1][2]:
                        merged_list.append(ranking[j])
                        count += 1
                        id_set.add(ranking[j][1])
                        j += 1
                    else:
                        break
                else:
                    j += 1

            if self.ranking[i][1] not in id_set: 
                if self.ranking[i] not in id_set and (count < self.args.topn or self.ranking[i][2] == merged_list[count-1][2]):
                    merged_list.append(self.ranking[i])  
                    count += 1
                    id_set.add(self.ranking[i][1])
                    i += 1
                else:
                    break
            else:
                i += 1

        while j < len(ranking):
            if ranking[j][1] not in id_set: 
                if ranking[j] not in id_set and (count < self.args.topn or ranking[j][2] == merged_list[count-1][2]):
                    merged_list.append(ranking[j])
                    count += 1
                    id_set.add(ranking[j][1])
                    j += 1
                else:
                    break
            else:
                j += 1

        self.ranking = merged_list

    def update_ranking(self, mol, max_tanimoto, ka_tag):
        index = 0
        if len(self.ranking) >= self.args.topn and max_tanimoto < self.ranking[len(self.ranking)-1][2]:
            pass
        else:    
            for top_mol in self.ranking:
                if max_tanimoto < top_mol[2]:
                    index = self.ranking.index(top_mol) + 1
                else:
                    break

            upper = self.ranking[:index]
            lower = self.ranking[index:]
            self.ranking = upper + [(oechem.OEMolToSmiles(mol), mol.GetTitle(), max_tanimoto, self.baitset[0], ka_tag)] + lower

            i = self.args.topn - 1
            while i < len(self.ranking) - 1:
                if self.ranking[i][2] != self.ranking[i + 1][2]:
                    self.ranking = self.ranking[:i + 1]

                    break
                else:
                    i += 1

class ParallelInsertKnownActives(ParallelComputeCube):
    """
    """

    classification = [["ParallelCompute"]]

    fptype = parameter.IntegerParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    data_input = ObjectInputPort('data_input')
    success = ObjectOutputPort('success')

    def process(self, data, port):
        
        self.act_list = data[0]
        self.baitset = data[1]
        self.ranking = data[2]
        self.fp_list = list()

        self.calculate_fp()
        self.insert_known_actives()

        self.success.emit((self.act_list, self.baitset, self.ranking))
        #self.success.emit(self.ranking)

    def calculate_fp(self):
        
        for mol in self.act_list:
            fp = oegraphsim.OEFingerPrint()
            oegraphsim.OEMakeFP(fp, mol, self.args.fptype)

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
        if len(self.ranking) >= self.args.topn and max_tanimoto < self.ranking[len(self.ranking)-1][2]:
            pass
        else:    
            for top_mol in self.ranking:
                if max_tanimoto < top_mol[2]:
                    index = self.ranking.index(top_mol) + 1
                else:
                    break

            upper = self.ranking[:index]
            lower = self.ranking[index:]
            self.ranking = upper + [(oechem.OEMolToSmiles(mol), mol.GetTitle(), max_tanimoto, self.baitset[0], ka_tag)] + lower

            i = self.args.topn - 1
            while i < len(self.ranking) - 1:
                if self.ranking[i][2] != self.ranking[i + 1][2]:
                    self.ranking = self.ranking[:i + 1]

                    break
                else:
                    i += 1

    def end(self):
        print('Parallel process ended')

class AccumulateRankings(ComputeCube):
    """
    A compute Cube that receives rankings and assemble them in a list
    """

    classification = [["Compute", "Accumulator"]]

    intake = ObjectInputPort('intake')
    success = ObjectOutputPort('success')

    def begin(self):
        self.ranking_list = list()

    def process(self, data, port):
        self.ranking_list.append(data[2])
        self.nb_ka = len(data[0])-len(data[1][1])

    def end(self):
        print('Accumulator ended')
        self.success.emit((self.ranking_list, self.nb_ka))

class AnalyseRankings(ComputeCube):
    """

    """

    classification = [["Compute", "Analysis"]]

    fptype = parameter.IntegerParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    intake = ObjectInputPort('intake')
    success = ObjectOutputPort('success')

    def process(self, data, port):
        self.ranking_list = data[0]
        self.nb_ka = data[1]
        print(str(self.nb_ka))
        self.results_avg = self.ranking_analysis()

        self.success.emit(self.results_avg)

    def ranking_analysis(self):
        results = pd.DataFrame()
        for i, ranking in enumerate(self.ranking_list):
            set_results = pd.DataFrame(columns = ['RR', 'HR', 'Set'])
            count = 0
            count_ka = 0
            for row, mol in enumerate(ranking):
                count += 1
                if mol[4] == 1:
                    count_ka += 1
                rr = 100 * count_ka/self.nb_ka
                hr = 100 * count_ka/count
                set_results.loc[row] = [rr, hr, i]
            results = pd.concat([results, set_results])
        
        results_avg = pd.DataFrame()

        fptypes = {102 : 'path', 104 : 'circular', 105 : 'tree'}
        FPType = fptypes[self.args.fptype]

        results_avg['Average RR ' + FPType] = results.groupby(results.index)['RR'].mean()
        results_avg['Average HR ' + FPType] = results.groupby(results.index)['HR'].mean()
        results_avg = results_avg.head(self.args.topn)

        return results_avg


class ParallelUpdateRanking(ParallelComputeCube):
    """
    A compute Cube that receives Molecules from the screening db with their Similarity value and rank them in the rnking list
    """

    classification = [["Compute", "Ranking"]]

    intake = ObjectInputPort('intake')
    success = ObjectOutputPort('success')

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    def begin(self):
        self.ranking= list()

    def process(self, data, port):
        index = 0
        if len(self.ranking) >= self.args.topn and data[2] < self.ranking[len(self.ranking)-1][2]:
            pass
        else:    
            for top_mol in self.ranking:
                if data[2] < top_mol[2]:
                    index = self.ranking.index(top_mol) + 1
                else:
                    break

            upper = self.ranking[:index]
            lower = self.ranking[index:]
            self.ranking = upper + [data] + lower

            i = self.args.topn - 1
            while i < len(self.ranking) - 1:
                if self.ranking[i][2] != self.ranking[i + 1][2]:
                    self.ranking = self.ranking[:i + 1]

                    break
                else:
                    i += 1

    def end(self):
        self.success.emit(self.ranking)

class ParallelCalculateFP(ParallelOEMolComputeCube):
    """
    A compute Cube that gets Molecules and assemble them in a list
    """

    classification = [["Compute", "Concat"]]

    fptype = parameter.IntegerParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")

    intake = MoleculeInputPort('intake')
    success = MoleculeOutputPort('success')


    def begin(self):
        pass

    def process(self, mol, port):
        fp = oegraphsim.OEFingerPrint()
        oegraphsim.OEMakeFP(fp, mol, self.args.fptype)

        time.sleep(random.random())

        self.success.emit(mol)

    def end(self):
        pass

class CreateTestData(ComputeCube):
    """
    Test Cube to delete
    """

    classification = [["Compute", "Test"]]

    fptype = parameter.IntegerParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")

    topn = parameter.IntegerParameter('topn', default=100,
                                    help_text="Number of top molecules returned in the rankinNumber of top molecules returned in the ranking")

    intake = ObjectInputPort('intake')
    success = ObjectOutputPort('success')

    def begin(self):

        fptypes = {102 : 'path', 104 : 'circular', 105 : 'tree'}
        FPType = fptypes[self.args.fptype]

        self.dataframe = pd.DataFrame(columns = ['Average RR ' + FPType, 'Average HR ' + FPType ])

    def process(self, data, port):
        self.set_id = data[0]
        self.baitset = data[1]
        for idx in self.baitset:
            self.dataframe.loc[len(self.dataframe)] = [idx, self.set_id]

    def end(self):
        self.success.emit(self.dataframe)



