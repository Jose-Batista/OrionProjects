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

import tempfile
from floe.api.orion import config_from_env, upload_file, stream_file
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

    intake = ObjectInputPort('intake')
    name = FileOutputParameter('name',
                               required=True,
                               description='The name of the output file')
    title = "Ranking Writer"
    classification = [["Output"]]

    def begin(self):
        pass

    def write(self, data, port):
        self.ranking_list = data[0]
        self.method = data[2]

        if self.method == 'Tree_FP':
            self.name_ext = 'FP_tree'
        elif self.method == 'FastROCS':
            self.name_ext = 'FR'

        self.in_orion = config_from_env() is not None
        if self.in_orion:
            self.stream = tempfile.NamedTemporaryFile()
        else:
            path = self.args.name + "ranking_" + self.name_ext + ".txt"
            self.stream = open(path, 'wb')
        for i, ranking in enumerate(self.ranking_list):
            text = "\n" + "Set n°" + str(ranking[0][2]) + "\n"
            text = text.encode("utf-8")
            self.stream.write(text)
            for mol in ranking:
                mol_data = str(mol[2]) + " " + mol[0] + " " + str(mol[1]) + " " + str(mol[3]) +  "\n"
                mol_data = mol_data.encode("utf-8")
                self.stream.write(mol_data)

        if self.in_orion:
            self.stream.flush()
            name = self.args.name + "ranking_" + self.name_ext + ".txt"
            resp = upload_file(name, self.stream.name)
            self.log.info("Created result file {} with ID {}".format(self.args.name, resp['id']))
        else:
            self.stream.close()

    def end(self):
        pass

class ResultsOutputCube(SinkCube):
    """
    A cube that outputs Results dataframe in a csv file
    """

    intake = ObjectInputPort('intake')
    name = FileOutputParameter('name',
                               required=True,
                               description='The name of the output file')
    title = "Results Writer"
    classification = [["Output"]]

    def begin(self):
        self.in_orion = config_from_env() is not None
        if self.in_orion:
            self.stream = tempfile.NamedTemporaryFile()

    def write(self, data, port):
        self.results_avg = data[0]
        self.method = data[1]

        if self.method == 'Tree_FP':
            self.name_ext = 'FP_tree'
        elif self.method == 'FastROCS':
            self.name_ext = 'FR'

        if self.in_orion:
            self.results_avg.to_csv(self.stream.name)
            self.stream.flush()
            name = self.args.name + "_results_" + self.name_ext + ".csv"
            resp = upload_file(name, self.stream.name)
            self.log.info("Created result file {} with ID {}".format(name, resp['id']))
        else:
            path = self.args.name + "results_" + self.name_ext + ".csv"
            self.results_avg.to_csv(path)

    def end(self):
        pass

class PlotResults(SinkCube):
    """

    """

    classification = [["Compute", "Plot"]]

    intake = ObjectInputPort('intake')
    name = FileOutputParameter('name',
                               required=True,
                               description='The name of the output file')

    def begin(self):
        self.results = pd.DataFrame()
        self.in_orion = config_from_env() is not None
        if self.in_orion:
            self.stream = tempfile.NamedTemporaryFile()

    def write(self, data, port):
        self.result = data[0]
        self.method = data[1]

        self.results = pd.concat([self.results, self.result], axis = 1)

#        if self.method == 'Tree_FP':
#            self.name_ext = 'FP_tree'
#        elif self.method == 'FastROCS':
#            self.name_ext = 'FR'
#
#        self.results_avg.plot(y = 'Average RR ' + self.name_ext, label = "Average RR" + self.name_ext)
#        plt.xlabel('Top Rank Molecules')
#        plt.ylabel('Rate (%)')
#        plt.legend( loc='best')
#        plt.title("Average RR Rates " + self.name_ext)
#        if self.in_orion:
#            plt.savefig(self.stream.name, format="svg")
#            self.stream.seek(0)
#            self.stream.flush()
#            name = self.args.name + "_Average_RR_plot_" + self.name_ext + ".svg"
#            resp = upload_file(name, self.stream.name)
#            self.log.info("Created result file {} with ID {}".format(name, resp['id']))
#        else:
#            path = self.args.name + "Average_RR_plot_" + self.name_ext + ".svg"
#            plt.savefig(path)
#
#        self.results_avg.plot(y = 'Average HR ' + self.name_ext, label = "Average HR" + self.name_ext)
#        plt.xlabel('Top Rank Molecules')
#        plt.ylabel('Rate (%)')
#        plt.legend( loc='best')
#        plt.title("Average HR Rates " + self.name_ext)
 
    def end(self):
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(9,8))
        plt.subplots_adjust(hspace=0.4)

        self.results.plot(ax=axes[0], y = [col for col in self.results.columns if 'RR' in col])
        axes[0].set_xlabel('Top Rank Molecules')
        axes[0].set_ylabel('Rate (%)')
        axes[0].set_title("Average RR Rates") 
        axes[0].legend([col.split(" ")[2] for col in self.results.columns if ('RR' in col)])
        
        self.results.plot(ax=axes[1], y = [col for col in self.results.columns if 'HR' in col])
        axes[1].set_xlabel('Top Rank Molecules')
        axes[1].set_ylabel('Rate (%)')
        axes[1].set_title("Average HR Rates") 
        axes[1].legend([col.split(" ")[2] for col in self.results.columns if ('HR' in col)])

        if self.in_orion:
            plt.savefig(self.stream.name, format="svg")
            self.stream.seek(0)
            self.stream.flush()
            name = self.args.name + "_Average_plot" + ".svg"
            resp = upload_file(name, self.stream.name)
            self.log.info("Created result file {} with ID {}".format(name, resp['id']))
        else:
            path = self.args.name + "Average_plot" + ".svg"
            plt.savefig(path)

class IndexOutputCube(SinkCube):
    """
    A cube that outputs an Index log
    """
    intake = ObjectInputPort('intake')
    name = FileOutputParameter('name',
                               required=True,
                               description='The name of the output file')
    title = "Index log Writer"
    tags = 'IndexOutputCube_tags' 
    classification = [["Output"]]

    def begin(self):
        self.in_orion = config_from_env() is not None
        if self.in_orion:
            self.stream = tempfile.NamedTemporaryFile()
        else:
            self.stream = open(self.args.name, 'wb')

    def write(self, data, port):
        self.set_id = data[0]
        self.baitset = data[1]
        text = 'Set n°' + str(self.set_id) + ': '
        text = text.encode('utf-8')
        self.stream.write(text)
        for idx in self.baitset:
            text = str(idx) + ' '
            text = text.encode('utf-8')
            self.stream.write(text)
        self.stream.write(b'\n')

    def end(self):
        if self.in_orion:
            self.stream.flush()
            resp = upload_file(self.args.name, self.stream.name)
            self.log.info("Created result file {} with ID {}".format(self.args.name, resp['id']))
        else:
            self.stream.close()

class TextOutputCube(SinkCube):
    """
    A cube that outputs text
    """
    intake = BinaryInputPort('intake')
    name = FileOutputParameter('name',
                               required=True,
                               description='The name of the output file')
    title = "Text File Writer"
    classification = [["Output"]]

    def begin(self):
        self.in_orion = config_from_env() is not None
        if self.in_orion:
            self.stream = tempfile.NamedTemporaryFile()
        else:
            self.stream = open(self.args.name, 'w')

    def write(self, data, port):
        data.decode('utf-8')
        self.stream.write(data)

    def end(self):
        if self.in_orion:
            self.stream.flush()
            resp = upload_file(self.args.name, self.stream.name)
            self.log.info("Created result file {} with ID {}".format(self.args.name, resp['id']))
        else:
            self.stream.close()

