import os

OUTDIR = 'results/scratch'
LOGDIR = 'logs'

[os.makedirs(x) for x in [OUTDIR, LOGDIR] if not os.path.exists(x)]

SIMS = list(range(1000-2000))
SIMS_PER_RUN = 100


rule all: 
	input: 'results/sim_aggregate.rdata'

rule sim_run:
	output: '{OUTDIR}/sim_{sim}.rdata'
	shell:
		'Rscript scripts/sim-run.r {wildcards.sim} {SIMS_PER_RUN} {output}'
rule sim_metrics:
	input: '{OUTDIR}/sim_{sim}.rdata'
	output: '{OUTDIR}/metrics_{sim}.rdata'
	shell:
		'Rscript scripts/sim-metrics.r {input} {output} {wildcards.sim}'

rule sim_aggregate:
	input: expand('{OUTDIR}/metrics_{sim}.rdata', OUTDIR=OUTDIR, sim=SIMS)
	output: 'results/sim_aggregate.rdata'
	shell:
		'Rscript scripts/sim-aggregate.r {input} {output}'

# rule analyse:
	
# rule learn:

# rule comparisons:

