import pickle

from openeye import oechem
from openeye import oeomega
from openeye import oemolprop

from vs_classes import VirtualScreeningData, ObjectInputPort

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
    intake = ObjectInputPort('intake')
    name = FileOutputParameter('name',
                               required=True,
                               description='The name of the output file')
    title = "File Writer"
    classification = [["Output"]]

    def begin(self):
        self.stream = open(self.args.name, 'w')

    def write(self, data, port):
        for i in range (len(data)):
            print(data[i].GetTitle())
            self.stream.write(data[i].GetTitle() + '\n')

    def end(self):
        self.stream.close()
