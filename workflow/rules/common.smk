import pandas as pd
import glob
import os

from snakemake.utils import validate

configfile: "../config/config.yaml"

validate(config, schema = "../schemas/config.schema.yaml")

units = pd.read_table(config["units"], dtype = str).set_index(
    ["sample"], drop = False
)
validate(units, schema = "../schemas/units.schema.yaml")

sample=units['sample']

outpath = config['outpath']

intervals_dir=config["interval"]

ref_version=config["ref_version"]

# Get the list of all interval files
interval_files = glob.glob(os.path.join(intervals_dir, "*.intervals.list"))

# Extract chromosome names by stripping the file path and extension
chromosomes = [os.path.basename(f).replace(".intervals.list", "") for f in interval_files]
cov=config['depth']

def get_unsort_bam(wildcards):
    """get fastq files of given sample"""
    bam = units.loc[(wildcards.sample), ["unsort_bam"]].dropna()
    return {"unsort_bam": bam.unsort_bam}
