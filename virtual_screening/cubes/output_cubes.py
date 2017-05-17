
from openeye import oechem
from openeye import oeomega
from openeye import oemolprop

from floe.api.parameter import (IntegerParameter, DataSetInputParameter, FileOutputParameter, FileInputParameter,
                                DataSetOutputParameter, BaseParameter, ParameterGroup,
                                DecimalParameter, StringParameter)

from floe.api.ports import (InputPort, OutputPort, Port, MoleculeInputPort,
                            MoleculeOutputPort, BinaryMoleculeInputPort, BinaryOutputPort,
                            MoleculeSerializerMixin, BinaryInputPort)

from floe.api import ( parameter, ParallelOEMolComputeCube, OEMolComputeCube, SinkCube)

class TextOutputCube(SinkCube):
    """
    A cube that outputs text
    """
    intake = BinaryInputPort('intake')
    name = FileOutputParameter('name',
                               required=True,
                               description='The name of the output file')
    title = "File Writer"
    classification = [["Output"]]

    def begin(self):
        self.stream = open(self.args.name, 'wb')

    def write(self, data, port):
        self.stream.write(data)

    def end(self):
        self.stream.close()
