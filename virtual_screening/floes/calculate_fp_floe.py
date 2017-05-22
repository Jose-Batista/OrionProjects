#!/usr/bin/env python

from cubes.input_cubes import IndexInputCube, OEMolTriggeredIStreamCube
from cubes.compute_cubes import CalculateFPCube, ParallelCalculateFP, ConcatMolList, GetSimValCube
from cubes.output_cubes import TextOutputCube
from floe.api import WorkFloe, CubeGroup
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('Read Text File')
job.description = """
Read an index text file and write the indices in an other file
"""
job.classification = [["Training", "Input/Output"]]

# Declare Cubes
act_reader = OEMolIStreamCube('act_reader')
act_reader.promote_parameter('data_in', promoted_name='act_db')
index_in = IndexInputCube('index_in')
index_in.promote_parameter('data_in', promoted_name='index_log')
#db_reader = OEMolTriggeredIStreamCube('db_reader')
#db_reader.promote_parameter('data_in', promoted_name='screen_db')

calculate_fp = CalculateFPCube('calculate_fp')
concat_act = ConcatMolList('concat_act')
get_sim_val = GetSimValCube('get_sim_val')
get_sim_val.promote_parameter('data_in', promoted_name='screen_db')

ofs = TextOutputCube('ofs')
ofs.promote_parameter('name', promoted_name='ofs')

# Create Cube group
group = CubeGroup(cubes=[calculate_fp, get_sim_val])

# Add Groups to Workfloe
job.add_group(group)

# Add Cubes to Floe
job.add_cubes(act_reader, index_in, calculate_fp, concat_act, get_sim_val, ofs)

# Connect ports
act_reader.success.connect(concat_act.intake)
concat_act.success.connect(calculate_fp.intake)
calculate_fp.success.connect(get_sim_val.fp_input)

index_in.success.connect(get_sim_val.baitset_input)

get_sim_val.success.connect(ofs.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()
