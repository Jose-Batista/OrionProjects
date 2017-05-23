#!/usr/bin/env python

from cubes.input_cubes import IndexInputCube, OEMolTriggeredIStreamCube
from cubes.compute_cubes import (CalculateFPCube, ParallelCalculateFP, ConcatMolList, ParallelRanking, ParallelUpdateRanking, 
                                PrepareRanking, ParallelInsertKnownActives, AccumulateRankings, AnalyseRankings)
from cubes.output_cubes import TextOutputCube, PlotResults
from floe.api import WorkFloe, CubeGroup
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('Read Text File')
job.description = """
Read an index text file and write the indices in an other file
"""
job.classification = [["Virtual Screening", "Create Ranking"]]

# Declare Cubes
act_reader = OEMolIStreamCube('act_reader')
act_reader.promote_parameter('data_in', promoted_name='act_db')
index_reader = IndexInputCube('index_reader')
index_reader.promote_parameter('data_in', promoted_name='index_log')

accu_act = ConcatMolList('accumulate actives')
#calc_fp = CalculateFPCube('calculate fingerprints')

prep_sim_calc = PrepareRanking('prepare similarity calculation')
calc_sim = ParallelRanking('calculate similarity value')
#calc_sim.promote_parameter('data_in', promoted_name='screen_db')

insert_known_actives = ParallelInsertKnownActives('insert known actives')
#update_ranking = ParallelUpdateRanking('update ranking')

accu_rankings = AccumulateRankings('accumulate rankings')
analyse_rankings = AnalyseRankings('analyse rankings')

plot_results = PlotResults('plot results')
plot_results.promote_parameter('name', promoted_name='plot_output')
ofs = TextOutputCube('ofs')
ofs.promote_parameter('name', promoted_name='ofs')

# Create Cube group
#group = CubeGroup(cubes=[prep_sim_calc, calc_sim])

# Add Groups to Workfloe
#job.add_group(group)

# Add Cubes to Floe
job.add_cubes(act_reader, index_reader, accu_act, prep_sim_calc, calc_sim, insert_known_actives, 
              accu_rankings, analyse_rankings, plot_results, ofs)

# Connect ports
act_reader.success.connect(accu_act.intake)
accu_act.success.connect(prep_sim_calc.act_input)
#calc_fp.success.connect(prep_sim_calc.fp_input)
index_reader.success.connect(prep_sim_calc.baitset_input)

prep_sim_calc.success.connect(calc_sim.data_input)
calc_sim.success.connect(insert_known_actives.data_input)
insert_known_actives.success.connect(accu_rankings.intake)

accu_rankings.success.connect(analyse_rankings.intake)
analyse_rankings.success.connect(plot_results.intake)
analyse_rankings.success.connect(ofs.intake)
#update_ranking.success.connect(ofs.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()

