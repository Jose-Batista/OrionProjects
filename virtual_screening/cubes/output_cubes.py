import pickle
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
        print(data)
        self.stream.write('test')

    def end(self):
        self.stream.close()

class PlotResults(SinkCube):
    """

    """

    classification = [["Compute", "Plot"]]

    fptype = parameter.IntegerParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")

    intake = ObjectInputPort('intake')
    name = FileOutputParameter('name',
                               required=True,
                               description='The name of the output file')

    def write(self, data, port):
        self.results_avg = data

        fptypes = {102 : 'path', 104 : 'circular', 105 : 'tree'}
        FPType = fptypes[self.args.fptype]

        self.results_avg.plot(y = 'Average RR ' + FPType, label = "Average RR" + FPType)
        plt.xlabel('Top Rank Molecules')
        plt.ylabel('Rate (%)')
        plt.legend( loc='best')
        plt.title("Average RR Rates " + FPType)
        path = self.args.name + "Average_ffp_RR_plot_" + FPType + ".svg"
        plt.savefig(path)

        self.results_avg.plot(y = 'Average HR ' + FPType, label = "Average HR" + FPType)
        plt.xlabel('Top Rank Molecules')
        plt.ylabel('Rate (%)')
        plt.legend( loc='best')
        plt.title("Average HR Rates FP" + FPType)
        path = self.args.name + "Average_ffp_HR_plot_" + FPType + ".svg"
        plt.savefig(path)
        
        plt.show()
