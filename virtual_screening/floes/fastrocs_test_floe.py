
#!/usr/bin/env python

from cubes.input_cubes import Test 
from cubes.compute_cubes import ParallelFastROCSRanking, AccumulateRankings
from cubes.output_cubes import TextRankingOutputCube
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
input_cube = Test('input')

request_cube = ParallelFastROCSRanking('request_cube')
request_cube.promote_parameter('url', promoted_name='url')
accu_cube = AccumulateRankings('accu')
accu_cube.promote_parameter('url', promoted_name='url')

output_cube = TextRankingOutputCube('results_output')
output_cube.promote_parameter('name', promoted_name='name')

job.add_cubes(input_cube, request_cube, accu_cube, output_cube)

input_cube.success.connect(request_cube.data_input)
request_cube.success.connect(accu_cube.rocsintake)
accu_cube.success.connect(output_cube.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()

