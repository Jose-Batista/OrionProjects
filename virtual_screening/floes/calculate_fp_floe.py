#!/usr/bin/env python

from cubes.input_cubes import IndexInputCube
from cubes.compute_cubes import CalculateFPCube
from cubes.output_cubes import TextOutputCube
from floe.api import WorkFloe 
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('Read Text File')
job.description = """
Read an index text file and write the indices in an other file
"""
job.classification = [["Training", "Input/Output"]]

# Declare Cubes
ifs = OEMolIStreamCube('ifs')
ifs.promote_parameter('data_in', promoted_name='ifs')
index_in = IndexInputCube('index_in')
index_in.promote_parameter('data_in', promoted_name='index_log')

calculate_fp = CalculateFPCube('calculate_fp')

ofs = TextOutputCube('ofs')
ofs.promote_parameter('name', promoted_name='ofs')

# Add Cubes to Floe
[job.add_cube(n) for n in [ifs, index_in, calculate_fp, ofs]]

# Connect ports
ifs.success.connect(calculate_fp.intake)
index_in.success.connect(calculate_fp.baitset_in)
calculate_fp.success.connect(ofs.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()
