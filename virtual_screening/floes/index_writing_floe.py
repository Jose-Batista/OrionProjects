
#!/usr/bin/env python

from cubes.compute_cubes import IndexGenerator, CreateTestData
from cubes.input_cubes import IndexInputCube
from cubes.output_cubes import IndexOutputCube, ResultsOutputCube, PlotResults
from floe.api import WorkFloe, CubeGroup
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('Test Outputs')
job.classification=[['Test']]
job.tags=[['funky cold medina'],]
job.title='test output'
job.description = """
Read an index text file and write outputs in several format
"""
input_cube = IndexInputCube('text_input')
input_cube.promote_parameter('data_in', promoted_name='Input')

create_test_data = CreateTestData('create test data')

#index_generator = IndexGenerator('index generator')
results_output = ResultsOutputCube('results output in csv')
results_output.promote_parameter('name', promoted_name='Output')
plot_results = PlotResults('plot results in svg')
plot_results.promote_parameter('name', promoted_name='Output')

index_output = IndexOutputCube('index output')
index_output.promote_parameter('name', promoted_name='Output')

job.add_cubes(input_cube, create_test_data, results_output, plot_results, index_output)
#job.add_cubes(input_cube, index_output)

input_cube.success.connect(index_output.intake)
input_cube.success.connect(create_test_data.intake)
create_test_data.success.connect(results_output.intake)
create_test_data.success.connect(plot_results.intake)
#index_generator.success.connect(index_output.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()
