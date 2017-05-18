
from openeye import oechem
from openeye import oeomega
from openeye import oemolprop

from floe.api.parameter import (IntegerParameter, DataSetInputParameter, FileOutputParameter, FileInputParameter,
                                DataSetOutputParameter, BaseParameter, ParameterGroup,
                                DecimalParameter, StringParameter)

from floe.api.ports import (InputPort, OutputPort, Port, MoleculeInputPort,
                            MoleculeOutputPort, BinaryMoleculeInputPort, BinaryOutputPort,
                            MoleculeSerializerMixin, BinaryInputPort)

from floe.constants import ADVANCED

from floe.api import (
    parameter, ParallelOEMolComputeCube, OEMolComputeCube, SourceCube,
    MoleculeOutputPort
)

class CalculateFPCube(ComputeCube):
    """
    A compute Cube that reads Molecules and calculate the fingerprint
    """

    classification = [["Compute", "Fingerprint"]]

    fptype = parameter.IntParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")
    success = OutputPort('success')
    limit = IntegerParameter('limit',
                             required=False,
                             level=ADVANCED,
                             description='Read up to N items from this cube')
    data_in = DataSetInputParameter('data_in',
                                    required=True,
                                    description='The index log to read from')

    def begin(self):

    def process(self, mol):
        fp = OEFingerprint()
        OEMakeFP(fp, mol, self.args.fptype)

        self.success.emit(mol.CreateCopy())

