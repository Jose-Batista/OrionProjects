import time
from floe.api import (
    OEMolComputeCube, SinkCube, BinaryInputPort, 
    InputPort, MoleculeInputPort, MoleculeOutputPort, ComputeCube
)

#from floe.api.orion import MultipartDatasetUploader, config_from_env

class TestCube(OEMolComputeCube):
#    molintake = MoleculeInputPort('molintake')
    txtintake = BinaryInputPort('txtintake')

    def begin(self):
        print("Hello World, TestCube is instantiated")

    def process(self, data, port):
        print( "Receiving something from port", port )
        print( "Data is", data )



