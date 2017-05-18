import pickle

from openeye import oechem
from openeye import oeomega
from openeye import oemolprop

from vs_classes import VirtualScreeningData, ObjectOutputPort

from floe.api import (parameter, ParallelOEMolComputeCube, OEMolComputeCube, SourceCube,
                            MoleculeOutputPort)

from floe.api.parameter import (IntegerParameter, DataSetInputParameter, FileOutputParameter, FileInputParameter,
                                DataSetOutputParameter, BaseParameter, ParameterGroup,
                                DecimalParameter, StringParameter)

from floe.api.ports import (InputPort, OutputPort, Port, MoleculeInputPort,
                            MoleculeOutputPort, BinaryMoleculeInputPort, BinaryOutputPort,
                            MoleculeSerializerMixin, BinaryInputPort)

class IndexInputCube(SourceCube):
    """
    An input cube that reads an index log and return the baitsets
    """

    classification = [["Input"]]

    success = ObjectOutputPort('success')
    limit = IntegerParameter('limit',
                             required=False,
                             description='Read up to N items from this cube')
    data_in = DataSetInputParameter('data_in',
                                    required=True,
                                    description='The index log to read from')

    data = VirtualScreeningData()

    def begin(self):
        self.stream = open(str(self.args.data_in), 'r')

    def __iter__(self):
        max_idx = self.args.limit
        if max_idx is not None:
            max_idx = int(max_idx)
        count = 0
        for line in self.stream:
            self.data.baitset.append(line)
            count += 1
            if max_idx is not None and count == max_idx:
                break
        yield self.data

