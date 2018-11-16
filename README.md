
```
conda env create -f environment.yml
conda activate r-simulations
snakemake -pr -j 100 \
--cluster-config bc4-cluster.json \
--cluster "sbatch \
--partition {cluster.partition} \
--nodes {cluster.nodes} \
--cpus-per-task {cluster.cpus-per-task} \
--time {cluster.time} \
--mem {cluster.mem} \
--output {cluster.output}"
```
