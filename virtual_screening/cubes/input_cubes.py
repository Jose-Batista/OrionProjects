
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

class IndexInputCube(SourceCube):
    """
    An input cube that reads an index log and return the baitsets
    """

    classification = [["Input"]]

    success = BinaryOutputPort('success')
    limit = IntegerParameter('limit',
                             required=False,
                             level=ADVANCED,
                             description='Read up to N items from this cube')
    data_in = DataSetInputParameter('data_in',
                                    required=True,
                                    description='The index log to read from')

    def begin(self):
        self.stream = open(str(self.args.data_in), 'rb')

    def __iter__(self):
        max_idx = self.args.limit
        if max_idx is not None:
            max_idx = int(max_idx)
        count = 0
        for line in self.stream:
            baitset = line.decode()
            baitset = baitset.split(" ")
            baitset = baitset[2:]
            print(baitset)
            for idx in baitset:
                idx = idx.encode()
                yield idx
            count += 1
            if max_idx is not None and count == max_idx:
                break


