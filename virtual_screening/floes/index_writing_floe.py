
#!/usr/bin/env python

from cubes.compute_cubes import IndexGenerator
from cubes.output_cubes import IndexOutputCube 
from floe.api import WorkFloe, CubeGroup
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('Demoxxx')
job.description = """
Read an index text file and write the indices in an other file
"""
index_generator = IndexGenerator('index generator')
index_output = IndexOutputCube('index output')
index_output.promote_parameter('data_in', promoted_name='index_file')

job.add_cubes(index_generator, index_output)

index_generator.success.connect(index_output.intake)
# If called from command line, run the floe
if __name__ == "__main__":
    job.run()
