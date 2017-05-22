
from openeye import oechem
from openeye import oegraphsim
from openeye import oeomega
from openeye import oemolprop

import time
import random

from vs_classes import VirtualScreeningData, ActiveList, ObjectInputPort, ObjectOutputPort

#from floe.api.parameter import (IntegerParameter, DataSetInputParameter, FileOutputParameter, FileInputParameter,
#                                DataSetOutputParameter, BaseParameter, ParameterGroup,
#                                DecimalParameter, StringParameter)

from floe.api.ports import (InputPort, OutputPort, Port, MoleculeInputPort,
                            MoleculeOutputPort, BinaryMoleculeInputPort, BinaryOutputPort,
                            MoleculeSerializerMixin, BinaryInputPort)

from floe.api import (parameter, ParallelOEMolComputeCube, OEMolComputeCube, SourceCube,
                      ComputeCube, ParallelComputeCube,
                      MoleculeOutputPort)


class ConcatMolList(ComputeCube):
    """
    A compute Cube that gets Molecules and assemble them in a list
    """

    classification = [["Compute", "Concat"]]

    intake = MoleculeInputPort('intake')
    success = ObjectOutputPort('success')

    def begin(self):
        self.mol_list = list()
        pass

    def process(self, mol, port):
        self.mol_list.append(mol)
        pass

    def end(self):
        self.success.emit(self.mol_list)
        

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
        pass

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
            pass

    def end(self):
  #      for fp in self.fp_list:
        self.success.emit(self.fp_list)
class PrepareSimCalc(ComputeCube):
    
    fp_input = ObjectInputPort('fp_input')
    baitset_input = ObjectInputPort('baitset_input')
    success = ObjectOutputPort('success')
    def begin(self):
        self.baitsets = list()
        self.fp_list = list()

    def process(self, data, port):
        if port is fp_input:
            self.fp_list = data

        if port is baitset_input:
            self.baitsets.append(data)

        if len(self.fp_list) > 0 and len(self.baitsets) > 0:
            self.success.emit((self.fp_list, self.baitsets.pop()))

class GetSimValCube(ComputeCube):
    """
    A compute Cube that receives a Molecule and a list of Fingerprints with a baitset of indices
    and returns the max Similarity value of the Molecule against the Fingerprints
    """

    classification = [["Compute", "Fingerprint", "Similarity"]]

    fptype = parameter.IntegerParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")

    data_in= parameter.DataSetInputParameter('data_in',
                                    required=True,
                                    title='Dataset to read from',
                                    description='The dataset to read from')
    download_format = parameter.StringParameter(
        'download_format',
        choices=('.oeb.gz', '.oeb', '.smi', '.pdb', '.mol2'),
        required=False,
        description='The stream format to be used for retrieving molecules from Orion',
        default=".oeb.gz")

    #intake = MoleculeInputPort('intake')
    act_data_input = ObjectInputPort('act_data_input')
    success = ObjectOutputPort('success')


    def begin(self):
#        self.max_tanimoto = 0
#        self.fp = None
        self.fp_list = None
        self.baitset = None
        pass

    def process(self, data, port):
       # if port == 'intake':
       #     print('mol intake')
       #     self.mol = 'test' #data.GetTitle()
       #     self.fp = oegraphsim.OEFingerPrint()
       #     oegraphsim.OEMakeFP(self.fp, data, self.args.fptype)

        self.fp_list = data[0]
        self.baitset = data[1]
        
        if self.fp_list is not None and self.baitset is not None:
            with oechem.oemolistream(str(self.args.data_in)) as ifs:
                for mol in ifs.GetOEMols():
                    max_tanimoto = 0
                    fp = oegraphsim.OEFingerPrint()
                    oegraphsim.OEMakeFP(fp, mol, self.args.fptype)
                    for idx in self.baitset:
                        act_fp = self.fp_list[idx]
                        tanimoto = oegraphsim.OETanimoto(fp, self.fp_list[idx])
                        if tanimoto > max_tanimoto:
                            max_tanimoto = tanimoto

                    self.success.emit((oechem.OEMolToSmiles(mol), mol.GetTitle(), max_tanimoto))
                    print('Molecule : Similarity %.3f' %(max_tanimoto))

        pass

    def end(self):
        pass

class UpdateRanking(ComputeCube):
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
            return ranking
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
