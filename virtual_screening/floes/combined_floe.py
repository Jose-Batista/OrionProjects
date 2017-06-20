
#!/usr/bin/env python

from cubes.input_cubes import IndexInputCube
from cubes.compute_cubes import (AccuMolList, PrepareRanking, ParallelFastFPRanking, ParallelFastROCSRanking,
                                ParallelFastFPInsertKA, ParallelInsertKARestfulROCS, AccumulateRankings, AnalyseRankings, IndexGenerator)
from cubes.output_cubes import TextRankingOutputCube, PlotResults, ResultsOutputCube
from floe.api import WorkFloe, CubeGroup
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('FastROCS VS Floe')
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
prep_ranking.promote_parameter('url', promoted_name='fastrocs_url')

create_FPranking = ParallelFastFPRanking('create_FPranking')
create_FPranking.promote_parameter('url', promoted_name='fastfp_url')
create_FPranking.promote_parameter('topn', promoted_name='topn')
create_ROCSranking = ParallelFastROCSRanking('create_ROCSranking')
create_ROCSranking.promote_parameter('url', promoted_name='fastrocs_url')
create_ROCSranking.promote_parameter('topn', promoted_name='topn')

insert_FPka = ParallelFastFPInsertKA('insert known actives in Fingerprint ranking')
insert_FPka.promote_parameter('url', promoted_name='fastfp_url')
insert_FPka.promote_parameter('topn', promoted_name='topn')
insert_ROCSka = ParallelInsertKARestfulROCS('insert known actives in FastROCS ranking')
insert_ROCSka.promote_parameter('url', promoted_name='fastrocs_url')
insert_ROCSka.promote_parameter('topn', promoted_name='topn')

accu_rankings = AccumulateRankings('accumulate rankings')
accu_rankings.promote_parameter('url', promoted_name='fastrocs_url')
analyse_rankings = AnalyseRankings('analyse rankings')
analyse_rankings.promote_parameter('fptype', promoted_name='fptype')
analyse_rankings.promote_parameter('topn', promoted_name='topn')

write_ranking = TextRankingOutputCube('write ranking')
write_ranking.promote_parameter('name', promoted_name='output_dir')
write_ranking.promote_parameter('fptype', promoted_name='fptype')
write_ranking.promote_parameter('method', promoted_name='method')
results_output = ResultsOutputCube('results output')
results_output.promote_parameter('name', promoted_name='output_dir')
results_output.promote_parameter('fptype', promoted_name='fptype')
plot_results = PlotResults('plot results')
plot_results.promote_parameter('name', promoted_name='output_dir')
plot_results.promote_parameter('fptype', promoted_name='fptype')


# Add Cubes to Floe
job.add_cubes(act_reader, index_generator, accu_act, prep_ranking, create_FPranking, create_ROCSranking, insert_FPka, 
              insert_ROCSka, accu_rankings, analyse_rankings, results_output, plot_results, write_ranking)

# Connect ports
act_reader.success.connect(accu_act.intake)
accu_act.success.connect(prep_ranking.act_input)
accu_act.success.connect(index_generator.intake)
index_generator.success.connect(prep_ranking.baitset_input)

prep_ranking.success.connect(create_FPranking.data_input)
prep_ranking.success.connect(create_ROCSranking.data_input)

create_FPranking.success.connect(insert_FPka.data_input)
create_ROCSranking.success.connect(insert_ROCSka.data_input)
insert_FPka.success.connect(accu_rankings.fpintake)
insert_ROCSka.success.connect(accu_rankings.rocsintake)

accu_rankings.success.connect(analyse_rankings.intake)

accu_rankings.success.connect(write_ranking.intake)
analyse_rankings.success.connect(results_output.intake)
analyse_rankings.success.connect(plot_results.intake)
#update_ranking.success.connect(write_ranking.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()

