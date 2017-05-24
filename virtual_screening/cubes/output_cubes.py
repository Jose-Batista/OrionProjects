import pickle
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from openeye import oechem
from openeye import oeomega
from openeye import oemolprop

from cubes.vs_classes import VirtualScreeningData, ObjectInputPort

from floe.api.parameter import (IntegerParameter, DataSetInputParameter, FileOutputParameter, FileInputParameter,
                                DataSetOutputParameter, BaseParameter, ParameterGroup,
                                DecimalParameter, StringParameter)

from floe.api.ports import (InputPort, OutputPort, Port, MoleculeInputPort,
                            MoleculeOutputPort, BinaryMoleculeInputPort, BinaryOutputPort,
                            MoleculeSerializerMixin, BinaryInputPort)

from floe.api import ( parameter, ParallelOEMolComputeCube, OEMolComputeCube, SinkCube)

class TextRankingOutputCube(SinkCube):
    """
    A cube that outputs text
    """

    fptype = parameter.IntegerParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")

    intake = ObjectInputPort('intake')
    name = FileOutputParameter('name',
                               required=True,
                               description='The name of the output file')
    title = "Ranking Writer"
    classification = [["Output"]]

    def begin(self):
        fptypes = {102 : 'path', 104 : 'circular', 105 : 'tree'}
        FPType = fptypes[self.args.fptype]

        path = self.args.name + "ranking_ffp_" + FPType + ".txt"
        self.stream = open(path, 'w')

    def write(self, data, port):
        self.ranking_list = data[0]

        for i, ranking in enumerate(self.ranking_list):
            self.stream.write("\n" + "Set n°" + str(ranking[0][3]) + "\n")
            for mol in ranking:
                mol_data = str(mol[3]) + " " + mol[1] + " " + str(mol[2]) + " " + str(mol[4]) +  "\n"
                self.stream.write(mol_data)

    def end(self):
        self.stream.close()

class ResultsOutputCube(SinkCube):
    """
    A cube that outputs Results dataframe in a csv file
    """

    fptype = parameter.IntegerParameter('fptype', default=105,
                                    help_text="Fingerprint type to use for the ranking")

    intake = ObjectInputPort('intake')
    name = FileOutputParameter('name',
                               required=True,
                               description='The name of the output file')
    title = "Results Writer"
    classification = [["Output"]]

    def begin(self):
        pass

    def write(self, data, port):
        self.results_avg = data

        fptypes = {102 : 'path', 104 : 'circular', 105 : 'tree'}
        FPType = fptypes[self.args.fptype]

        path = self.args.name + "results_ffp_" + FPType + ".csv"
        self.results_avg.to_csv(path)

    def end(self):
        pass

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

class IndexOutputCube(SinkCube):
    """
    A cube that outputs an Index log
    """
    intake = ObjectInputPort('intake')
    name = FileOutputParameter('name',
                               required=True,
                               description='The name of the output file')
    title = "Index log Writer"
    classification = [["Output"]]

    def begin(self):
        self.in_orion = config_from_env() is not None
        if self.in_orion:
            self.stream = tempfile.NamedTemporaryFile()
        else:
            self.stream = open(self.args.name, 'w')

    def write(self, data, port):
        self.set_id = data[0]
        self.baitset = data[1]
        self.stream.write('Set n°' + str(self.set_id) + ': ')
        for idx in self.baitset:
            self.stream.write(str(idx))
        self.stream.write('\n')

    def end(self):
        if self.in_orion:
            self.stream.flush()
            resp = upload_file(self.args.name, self.stream.name)
            self.log.info("Created result file {} with ID {}".format(self.args.name, resp['id']))
        else:
            self.stream.close()

