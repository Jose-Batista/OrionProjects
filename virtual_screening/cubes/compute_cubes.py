
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
    success = BinaryOutputPort('success')


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
        for fp in self.fp_list:
            self.success.emit(fp)

class GetSimValCube(ComputeCube):
    """
    A compute Cube that receives a Molecule and a list of Fingerprints with a baitset of indices
    and returns the max Similarity value of the Molecule against the Fingerprints
    """

    classification = [["Compute", "Fingerprint", "Similarity"]]

    fptype = parameter.IntegerParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")

    #intake = parameter.DataSetInputParameter('intake',
    #                                required=True,
    #                                title='Dataset to read from',
    #                                description='The dataset to read from')
    #download_format = StringParameter(
    #    'download_format',
    #    choices=('.oeb.gz', '.oeb', '.smi', '.pdb', '.mol2'),
    #    required=False,
    #    description='The stream format to be used for retrieving molecules from Orion',
    #    default=".oeb.gz")

    intake = MoleculeInputPort('intake')
    fp_input = BinaryInputPort('fp_input')
    baitset_input = ObjectInputPort('baitset_input')
    success = ObjectOutputPort('success')


    def begin(self):
        self.max_tanimoto = 0
        self.fp = None
        self.fp_list = None
        self.baitset = None
        pass

    def process(self, data, port):
        if port == 'intake':
            print('mol intake')
            self.mol = 'test' #data.GetTitle()
            self.fp = oegraphsim.OEFingerPrint()
            oegraphsim.OEMakeFP(self.fp, data, self.args.fptype)

        if port == 'fp_input':
            print('fp_list input')
            self.fp_list.append(data)

        if port == 'baitset_input':
            print('baitset input')
            self.baitset = data
        
        if self.fp is not None and self.fp_list is not None and self.baitset is not None:
            for idx in self.baitset:
                act_fp = self.fp_list[idx]
                tanimoto = oegraphsim.OETanimoto(self.fp, self.fp_list[idx])
                if tanimoto > self.max_tanimoto:
                    self.max_tanimoto = tanimoto

            self.success.emit(self.max_tanimoto)
            print('Molecule : Similarity %.3f' %(self.max_tanimoto))

        pass

    def end(self):
        pass

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
