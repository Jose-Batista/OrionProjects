
#!/usr/bin/env python

from cubes.compute_cubes import ShapeDatabaseClient
from cubes.input_cubes import IndexInputCube
from cubes.output_cubes import IndexOutputCube, ResultsOutputCube, PlotResults
from floe.api import WorkFloe, CubeGroup
from floe.api import OEMolIStreamCube
from floe.api import OEMolOStreamCube, FileOutputCube

# Declare Floe, add metadata for UI
job = WorkFloe('FastROCS Test Floe')
job.classification=[['Test']]
job.tags=[['yippee ki yay mf']]
job.title='test FastROCS Server'
job.description = """
Read a molecule query and return the FastROCS Server Results
"""
input_cube = OEMolIStreamCube('mol_input')

request_cube = ShapeDatabaseClient('request cube')

output_cube = FileOutputCube('results output')
output_cube.promote_parameter('name', promoted_name='output')

job.add_cubes(input_cube, request_cube, output_cube)

input_cube.success.connect(request_cube.intake)
request_cube.success.connect(output_cube.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()

