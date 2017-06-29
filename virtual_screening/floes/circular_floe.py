
#!/usr/bin/env python

from cubes.input_cubes import IndexInputCube
from cubes.compute_cubes import (AccuMolList, PrepareRanking, ParallelTreeFPRanking, ParallelPathFPRanking, ParallelCircularFPRanking, ParallelFastROCSRanking,
                                ParallelTreeFPInsertKA, ParallelPathFPInsertKA, ParallelCircularFPInsertKA, ParallelInsertKARestfulROCS, AccumulateRankingsCircular, AnalyseRankings, IndexGenerator)
from cubes.output_cubes import WriteRanking, PlotResults, ResultsOutputCube
from floe.api import WorkFloe, CubeGroup, ParallelCubeGroup
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('Circular Floe')
job.description = """
Virtual Screening of a database against actives with the FastROCS method
"""
job.classification = [["Virtual Screening", "Create Ranking"]]

# Declare Cubes
act_reader = OEMolIStreamCube('act_reader')
act_reader.promote_parameter('data_in', promoted_name='act_db')
#index_reader = IndexInputCube('index_reader')
#index_reader.promote_parameter('data_in', promoted_name='index_log')

index_generator = IndexGenerator('index generator')
accu_act = AccuMolList('accumulate actives')

prep_ranking = PrepareRanking('prepare similarity calculation')

create_circularFPranking = ParallelCircularFPRanking('create_circularFPranking')
create_circularFPranking.promote_parameter('url', promoted_name='fastfp_url')
create_circularFPranking.promote_parameter('topn', promoted_name='topn')

insert_circularFPka = ParallelCircularFPInsertKA('insert known actives in Circular Fingerprint ranking')
insert_circularFPka.promote_parameter('topn', promoted_name='topn')

accu_rankings = AccumulateRankingsCircular('accumulate rankings')
analyse_rankings = AnalyseRankings('analyse rankings')
analyse_rankings.promote_parameter('topn', promoted_name='topn')

write_ranking = WriteRanking('write ranking')
write_ranking.promote_parameter('name', promoted_name='output_dir')
results_output = ResultsOutputCube('results output')
results_output.promote_parameter('name', promoted_name='output_dir')
plot_results = PlotResults('plot results')
plot_results.promote_parameter('name', promoted_name='output_dir')

# Create Cube group
#tree_group = ParallelCubeGroup(cubes=[create_treeFPranking, insert_treeFPka])
#path_group = ParallelCubeGroup(cubes=[create_pathFPranking, insert_pathFPka])
#circular_group = ParallelCubeGroup(cubes=[create_circularFPranking, insert_circularFPka])
#rocs_group = ParallelCubeGroup(cubes=[create_ROCSranking, insert_ROCSka ])

# Add Groups to Workfloe
#for group in [tree_group, path_group, circular_group, rocs_group]:
#    job.add_group(group)
#
## Add Cubes to Floe
job.add_cubes(act_reader, index_generator, accu_act, prep_ranking, create_circularFPranking, 
              insert_circularFPka, accu_rankings, analyse_rankings, results_output, plot_results, write_ranking)

# Connect ports
act_reader.success.connect(accu_act.intake)
accu_act.success.connect(prep_ranking.act_input)
accu_act.success.connect(index_generator.intake)
index_generator.success.connect(prep_ranking.baitset_input)
#index_reader.success.connect(prep_ranking.baitset_input)

prep_ranking.success.connect(create_circularFPranking.data_input)

create_circularFPranking.success.connect(insert_circularFPka.data_input)
insert_circularFPka.success.connect(accu_rankings.circular_fpintake)

accu_rankings.success.connect(analyse_rankings.intake)

accu_rankings.success.connect(write_ranking.intake)
analyse_rankings.success.connect(results_output.intake)
analyse_rankings.success.connect(plot_results.intake)
#update_ranking.success.connect(write_ranking.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()


