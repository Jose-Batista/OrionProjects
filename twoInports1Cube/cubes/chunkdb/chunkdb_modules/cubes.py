import time
from floe.api import (
    OEMolComputeCube, DecimalParameter, JsonInputPort, JsonOutputPort,
    MoleculeInputPort, MoleculeOutputPort, ComputeCube
)


# For parallel, import and inherit from ParallelOEMolComputeCube
class SleepyCube(OEMolComputeCube):
    """
    docstring for your cube, only visible to developers

    This cube just sleeps
    """
    title = "Cube Title for view in UI"
    description = """
    *Longform Description*
    Shown in the UI for floe builders. Should explain cube and params
    """
    classification = [
        ["Testing", "Sleep"],
        ["Testing", "Example"],
    ]
    tags = [tag for lists in classification for tag in lists]

    amount = DecimalParameter("amount",
            title="Amount to Sleep",
            description="The amount to sleep (per mol)",
            default=1.0,
            )

    """
    For each molecule received the process function will be called
    """
    def process(self, mol, port):
        time.sleep(float(self.args.amount))
        self.emit(mol)


class AccumulatorCube(OEMolComputeCube):
    """
    Cube where you need to accumulate results so as to
    perform some sort of computation on a group of molecules
    or some other data
    """

    title = "Example Accumulator Cube"
    description = """
    An example of a cube that acts upon a number of Molecules at once

    Cannot Parallelize accumulators
    """
    classification = [
        ["Testing", "Accumulator"],
        ["Testing", "Example"],
    ]
    tags = [tag for lists in classification for tag in lists]

    """
    Initial setup of cube

    Create an array to store data that the cube receives
    """
    def begin(self):
        self.accumulated_data = []

    """
    For each OEMol received store the data in the array
    created in the begin() function
    """
    def process(self, data, port):
        self.accumulated_data.append(data)

    """
    Act upon the accumulated molecules

    Just loop through molecules and print title before emitting
    """
    def end(self):
        for mol in self.accumulated_data:
            print(mol.GetTitle())
            self.success.emit(mol)


class ExtraPortsCube(ComputeCube):
    """
    Cube where you define different types of ports

    Inherits from ComputeCube which has no default Ports like
    OEMolComputeCube
    """

    title = "Example Extra Ports Cube"
    description = """
    An Example of a cube that has both Molecule and Json ports defined
    on them and how to handle the different inputs in the process function
    """

    classification = [
        ["Testing", "Ports"],
        ["Testing", "Example"],
    ]
    tags = [tag for lists in classification for tag in lists]

    mol_in = MoleculeInputPort("mol_in")
    mol_out = MoleculeOutputPort("mol_out")

    json_in = JsonInputPort("json_in")
    json_out = JsonOutputPort("json_out")

    """
    For each OEMol received store the data in the array
    created in the begin() function
    """
    def process(self, data, port):
        if port == "mol_in":
            self.process_mol(data)
        elif port == "json_in":
            self.process_json(data)
        else:
            self.log.warn("Data received on unexpected Port")
            return

    def process_mol(self, mol):
        print(mol.GetTitle())
        self.mol_out.emit(mol)

    def process_json(self, data):
        for key in data:
            print(key)
        self.json_out.emit(data)
