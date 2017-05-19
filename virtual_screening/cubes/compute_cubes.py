
from openeye import oechem
from openeye import oegraphsim
from openeye import oeomega
from openeye import oemolprop

import time
import random

from vs_classes import VirtualScreeningData, ObjectInputPort, ObjectOutputPort

#from floe.api.parameter import (IntegerParameter, DataSetInputParameter, FileOutputParameter, FileInputParameter,
#                                DataSetOutputParameter, BaseParameter, ParameterGroup,
#                                DecimalParameter, StringParameter)

from floe.api.ports import (InputPort, OutputPort, Port, MoleculeInputPort,
                            MoleculeOutputPort, BinaryMoleculeInputPort, BinaryOutputPort,
                            MoleculeSerializerMixin, BinaryInputPort)

from floe.api import (parameter, ParallelOEMolComputeCube, OEMolComputeCube, SourceCube,
                      ComputeCube, ParallelComputeCube,
                      MoleculeOutputPort)


class ConcatActiveList(ComputeCube):
    """
    A compute Cube that gets Molecules and assemble them in a list
    """

    classification = [["Compute", "Concat"]]

    intake = MoleculeInputPort('intake')
    success = ObjectOutputPort('success')

    def begin(self):
        self.act_list = list()
        pass

    def process(self, mol, port):
        self.act_list.append(mol)
        pass

    def end(self):
        self.success.emit(self.act_list)
        

class CalculateFPCube(ComputeCube):
    """
    A compute Cube that reads Molecules and calculate the fingerprint
    """

    classification = [["Compute", "Fingerprint"]]

    fptype = parameter.IntegerParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")

    intake = ObjectInputPort('intake')
    baitset_in = ObjectInputPort('baitset_in')
    success = OutputPort('success')


    def begin(self):
        pass

    def process(self, mol, port):
        fp = OEFingerPrint()
        OEMakeFP(fp, mol, self.args.fptype)

        self.success.emit(mol.CreateCopy())


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
