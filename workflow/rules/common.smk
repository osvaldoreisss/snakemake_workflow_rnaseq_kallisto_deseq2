from snakemake.utils import validate
import pandas as pd
import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
#validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep=",").set_index(["sample","run"], drop=False)
#validate(samples, schema="../schemas/samples.schema.yaml")

threads = config['threads']

samples_file = config['samples']

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")

HTTP = HTTPRemoteProvider()

def get_fastq(wildcards):
    fastqs = samples.loc[(wildcards.sample, int(wildcards.run)), ["fq1", "fq2"]].dropna()
    return {'fq1': HTTP.remote(f"{fastqs.fq1}", keep_local=True), 'fq2': HTTP.remote(f"{fastqs.fq2}", keep_local=True)}
    