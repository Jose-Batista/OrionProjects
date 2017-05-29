
#!/usr/bin/env python

from cubes.compute_cubes import IndexGenerator
from cubes.input_cubes import IndexInputCube
from cubes.output_cubes import IndexOutputCube 
from floe.api import WorkFloe, CubeGroup
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('Demoxxx')
job.classification=[['blah']]
job.tags=[['funky cold medina'],]
job.title='mytitle'
job.description = """
Read an index text file and write the indices in an other file
"""
input_cube = IndexInputCube('text_input')
#input_cube.promote_parameter('data_in', promoted_name='index_input')
#index_generator = IndexGenerator('index generator')
index_output = IndexOutputCube('index output')
index_output.promote_parameter('name', promoted_name='index_output')

job.add_cubes(input_cube, index_output)

input_cube.success.connect(index_output.intake)
#index_generator.success.connect(index_output.intake)
# If called from command line, run the floe
if __name__ == "__main__":
    job.run()
