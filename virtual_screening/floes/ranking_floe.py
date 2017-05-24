#!/usr/bin/env python

from cubes.input_cubes import IndexInputCube, OEMolTriggeredIStreamCube
from cubes.compute_cubes import (CalculateFPCube, ParallelCalculateFP, AccuMolList, ParallelRanking, ParallelUpdateRanking, 
                                PrepareRanking, ParallelInsertKnownActives, AccumulateRankings, AnalyseRankings)
from cubes.output_cubes import TextRankingOutputCube, PlotResults, ResultsOutputCube
from floe.api import WorkFloe, CubeGroup
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('VS Floe')
job.description = """
Read an index text file and write the indices in an other file
"""
job.classification = [["Virtual Screening", "Create Ranking"]]

# Declare Cubes
act_reader = OEMolIStreamCube('act_reader')
act_reader.promote_parameter('data_in', promoted_name='act_db')
index_reader = IndexInputCube('index_reader')
index_reader.promote_parameter('data_in', promoted_name='index_log')

accu_act = AccuMolList('accumulate actives')
#calc_fp = CalculateFPCube('calculate fingerprints')

prep_sim_calc = PrepareRanking('prepare similarity calculation')
calc_sim = ParallelRanking('calculate similarity value')
calc_sim.promote_parameter('fptype', promoted_name='fptype')
calc_sim.promote_parameter('topn', promoted_name='topn')
#calc_sim.promote_parameter('data_in', promoted_name='screen_db')

insert_known_actives = ParallelInsertKnownActives('insert known actives')
insert_known_actives.promote_parameter('fptype', promoted_name='fptype')
insert_known_actives.promote_parameter('topn', promoted_name='topn')
#update_ranking = ParallelUpdateRanking('update ranking')

accu_rankings = AccumulateRankings('accumulate rankings')
analyse_rankings = AnalyseRankings('analyse rankings')
analyse_rankings.promote_parameter('fptype', promoted_name='fptype')
analyse_rankings.promote_parameter('topn', promoted_name='topn')

write_ranking = TextRankingOutputCube('write ranking')
write_ranking.promote_parameter('name', promoted_name='output_dir')
write_ranking.promote_parameter('fptype', promoted_name='fptype')
results_output = ResultsOutputCube('results output')
results_output.promote_parameter('name', promoted_name='output_dir')
results_output.promote_parameter('fptype', promoted_name='fptype')
plot_results = PlotResults('plot results')
plot_results.promote_parameter('name', promoted_name='output_dir')
plot_results.promote_parameter('fptype', promoted_name='fptype')

# Create Cube group
#group = CubeGroup(cubes=[prep_sim_calc, calc_sim])

# Add Groups to Workfloe
#job.add_group(group)

# Add Cubes to Floe
job.add_cubes(act_reader, index_reader, accu_act, prep_sim_calc, calc_sim, insert_known_actives, 
              accu_rankings, analyse_rankings, results_output, plot_results, write_ranking)

# Connect ports
act_reader.success.connect(accu_act.intake)
accu_act.success.connect(prep_sim_calc.act_input)
#calc_fp.success.connect(prep_sim_calc.fp_input)
index_reader.success.connect(prep_sim_calc.baitset_input)

prep_sim_calc.success.connect(calc_sim.data_input)
calc_sim.success.connect(insert_known_actives.data_input)
insert_known_actives.success.connect(accu_rankings.intake)

accu_rankings.success.connect(analyse_rankings.intake)

accu_rankings.success.connect(write_ranking.intake)
analyse_rankings.success.connect(results_output.intake)
analyse_rankings.success.connect(plot_results.intake)
#update_ranking.success.connect(write_ranking.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()

