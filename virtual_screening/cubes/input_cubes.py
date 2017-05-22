import pickle

from openeye import oechem
from openeye import oeomega
from openeye import oemolprop

from vs_classes import VirtualScreeningData, ObjectOutputPort, ObjectInputPort

from floe.api import (parameter, ParallelOEMolComputeCube, OEMolComputeCube, SourceCube, ComputeCube,
                      MoleculeOutputPort)

#from floe.api.parameter import (IntegerParameter, DataSetInputParameter, FileOutputParameter, FileInputParameter,
#                                DataSetOutputParameter, BaseParameter, ParameterGroup,
#                                DecimalParameter, StringParameter)

from floe.api.ports import (InputPort, OutputPort, Port, MoleculeInputPort,
                            MoleculeOutputPort, BinaryMoleculeInputPort, BinaryOutputPort,
                            MoleculeSerializerMixin, BinaryInputPort)

class IndexInputCube(SourceCube):
    """
    An input cube that reads an index log and return the baitsets
    """

    classification = [["Input"]]

    success = ObjectOutputPort('success')

    limit = parameter.IntegerParameter('limit',
                             required=False,
                             description='Read up to N items from this cube')

    data_in = parameter.DataSetInputParameter('data_in',
                                    required=True,
                                    description='The index log to read from')


    def begin(self):
        #self.baitset = list()
        self.stream = open(str(self.args.data_in), 'r')

    def __iter__(self):
        max_idx = self.args.limit
        if max_idx is not None:
            max_idx = int(max_idx)
        count = 0

        for line in self.stream:
            self.baitset = line
            self.baitset = self.baitset.split(" ")
            self.baitset = self.baitset[2:-1]
            for i, idx in enumerate(self.baitset):
                self.baitset[i] = int(idx)
            count += 1
            if max_idx is not None and count == max_idx:
                break
        yield self.baitset


class OEMolTriggeredIStreamCube(ComputeCube):
    """
    A source cube that uses oechem to read molecules
    """
    classification = [["Input"]]
    success = MoleculeOutputPort('success')

    title = "Dataset Reader"

    limit = parameter.IntegerParameter('limit',
                             required=False,
                             description='Read up to N items from this cube')
    fp_input = ObjectInputPort('fp_input')
    data_in = parameter.DataSetInputParameter('data_in',
                                    required=True,
                                    title='Dataset to read from',
                                    description='The dataset to read from')
    download_format = parameter.StringParameter(
        'download_format',
        choices=('.oeb.gz', '.oeb', '.smi', '.pdb', '.mol2'),
        required=False,
        description='The stream format to be used for retrieving molecules from Orion',
        default=".oeb.gz")

    received_act = False

    def process(self, data, port):
        #print(data,port)
        if port is 'fp_input':
            print('Curry wurst')
            self.received_act = True
            max_idx = self.args.limit
            if max_idx is not None:
                max_idx = int(max_idx)
            count = 0
            with oechem.oemolistream(str(self.args.data_in)) as ifs:
                for mol in ifs.GetOEMols():
                    self.success.emit(mol)
                    count += 1
                    if max_idx is not None and count == max_idx:
                        break
            #else:
            #    stream = StreamingDataset(self.args.data_in,
            #                              input_format=self.args.download_format)
            #    for mol in stream:
            #        count += 1
            #        self.success.emit(mol)
            #        if max_idx is not None and count == max_idx:
            #            break

#    def __iter__(self):
#        if self.received_act:
#            max_idx = self.args.limit
#            if max_idx is not None:
#                max_idx = int(max_idx)
#            count = 0
#            self.config = config_from_env()
#            in_orion = self.config is not None
#            if not in_orion:
#                with oemolistream(str(self.args.data_in)) as ifs:
#                    for mol in ifs.GetOEMols():
#                        print(mol.GetTitle())
#                        yield mol
#                        count += 1
#                        if max_idx is not None and count == max_idx:
#                            break
#            else:
#                stream = StreamingDataset(self.args.data_in,
#                                          input_format=self.args.download_format)
#                for mol in stream:
#                    yield mol
#                    count += 1
#                    if max_idx is not None and count == max_idx:
#                        break
#        else:
#            print('Waiting...')
#
