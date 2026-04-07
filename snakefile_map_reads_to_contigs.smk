import pandas as pd
import collections
import os
from pathlib import Path
import sys
import itertools
from snakemake.io import glob_wildcards

# Set the directory the snakefile exists in. This makes us able to call the pipeline with the relevant src files from other directories.
THIS_FILE_DIR = Path(workflow.basedir)

# If the configfile is not set explicit fall back to the default
CONFIG_PATH = THIS_FILE_DIR / "config/config.yaml"
# TODO should go to configfile
if CONFIG_PATH.exists():
    configfile: CONFIG_PATH

# Define the src directory for the files used in the snakemake workflow
SRC_DIR = THIS_FILE_DIR / "files_used_in_snakemake_workflow"

# Get the output_directory defined by the user or fallback to current directory, which is the default way snakemake handles output directories
OUTDIR = Path("") if config.get("output_directory") is None else Path(config.get("output_directory"))
if not OUTDIR.exists():
    OUTDIR.mkdir()

#### Setting parameters from the config file ####
##  For a more throughout description of what the different config options mean see the /config/config.yaml file

# Default resources used
default_walltime = config.get("default_walltime", "48:00:00")
default_threads = config.get("default_threads", 16)
default_mem_gb = config.get("default_mem_gb", 50)

# Minimum contig length used
MIN_CONTIG_LEN = int(config.get("min_contig_len", 2000))

# Binning parameters
PLAMB_PARAMS = config.get("plamb_params", ' -o C ')
PLAMB_PRELOAD = config.get("plamb_preload", "")

# Create opton for user to pass in vamb arguments from CLI. Usefull for easy example on small dataset.
vamb_arguments = config.get("vamb_arguments", None)
if vamb_arguments is not None:
    PLAMB_PARAMS = vamb_arguments

# Other options
CUDA = True if config.get("cuda") ==  "True" else False

## ----------- ##

# Assert that input files are actually passed to snakemake
if config.get("read_file") == None and config.get("read_assembly_dir") == None and config.get("should_install_genomad") == None:
    print("ERROR: read_file or read_assembly_dir not passed to snakemake as config. Define either. Eg. snakemake <arguments> --config read_file=<read file>. If in doubt refer to the README.md file")
    sys.exit()


# Set default values for dictonaries containg information about the input information
# The way snakemake parses snakefiles means we have to define them even though they will always be present
sample_id = dict()
sample_id_path= dict()
sample_id_path_assembly = dict()

# If the read_file is defined the pipeline will also run SPades and assemble the reads
if config.get("read_file") != None:
    df = pd.read_csv(config["read_file"], sep=r"\s+", comment="#")
    sample_id = collections.defaultdict(list)
    sample_id_path = collections.defaultdict(dict)
    for id, (read1, read2) in enumerate(zip(df.read1, df.read2)):
        id = f"sample{str(id)}"
        # Earlier version of the pipeline could handle passing several samples at the same time which would be processed separatly.
        # For easier user input this was removed. Therefore the "sample" is always set to the same.
        # By parsing different sample names from the input this can be implemented again - this is nice eg. for benchmarking
        sample = "intermidiate_files"
        sample_id[sample].append(id)
        sample_id_path[sample][id] = [read1, read2]

# If read_assembly dir is defined the pipeline will run user defined SPades output files
if config.get("read_assembly_dir") != None:
    df = pd.read_csv(config["read_assembly_dir"], sep=r"\s+", comment="#")
    sample_id = collections.defaultdict(list)
    sample_id_path = collections.defaultdict(dict)
    sample_id_path_assembly = collections.defaultdict(dict)
    for id, (read1, read2, assembly) in enumerate(zip( df.read1, df.read2, df.assembly_dir)):
        #id = f"sample{str(id)}"
        id = f"{str(id)}"
        sample = "intermidiate_files"
        sample_id[sample].append(id)
        sample_id_path[sample][id] = [read1, read2]
        sample_id_path_assembly[sample][id] = [assembly]

    # Setting the output paths for the user defined SPades files
    contigs =  lambda wildcards: Path(sample_id_path_assembly["intermidiate_files"][wildcards.id][0]) / "contigs.fasta"
    assembly_graph  =  lambda wildcards: Path(sample_id_path_assembly["intermidiate_files"][wildcards.id][0]) / "assembly_graph_after_simplification.gfa"
    contigs_paths  =  lambda wildcards: Path(sample_id_path_assembly["intermidiate_files"][wildcards.id][0]) / "contigs.paths"

read_fw = lambda wildcards: sample_id_path["intermidiate_files"][wildcards.id][0]
read_rv =  lambda wildcards: sample_id_path["intermidiate_files"][wildcards.id][1]



# Functions to get the config-defined threads/walltime/mem_gb for a rule and if not defined the default
threads_fn = lambda rulename: config.get(rulename, {"threads": default_threads}).get("threads", default_threads)
walltime_fn  = lambda rulename: config.get(rulename, {"walltime": default_walltime}).get("walltime", default_walltime)
mem_gb_fn  = lambda rulename: config.get(rulename, {"mem_gb": default_mem_gb}).get("mem_gb", default_mem_gb)

try:
    os.makedirs(os.path.join(OUTDIR,'log'), exist_ok=True)
except FileExistsError:
    pass

rulename = "all"
rule all:
    input:
        finished = os.path.join(OUTDIR,"intermidiate_files",'rule_completed_checks/run_VAE.finished')

# Run strobealign to get the abundances
rulename = "Strobealign_bam_default"
rule Strobealign_bam_default:
        input:
            fw = read_fw,
            rv = read_rv,
            contigs = "/home/projects/cu_10108/people/paupie/plasmaag_bangladesh_all_samples_wo_birth_1_2_2_2_based_on_louvain/results/candidate_plasmids_and_phage_contigs.fna"
        output:
            OUTDIR / "intermidiate_files/assembly_mapping_output/mapped/{id}.bam"
        threads: threads_fn(rulename)
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        benchmark: config.get("benchmark", f"{str(OUTDIR)}/benchmark/") + "intermidiate_files_{id}_" + rulename
        log: 
            log=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_{id}_" + rulename,
            e=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_{id}_" + rulename+"_err",
            o=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_{id}_" + rulename+"_out"
        conda: THIS_FILE_DIR / "envs/strobe_env.yaml"
        shell:
            """
            strobealign -t {threads} {input.contigs} {input.fw} {input.rv} > {output} 2> {log.log}
            """

# Sort the bam files and index them
rulename="sort"
rule sort:
    input:
        OUTDIR / "intermidiate_files/assembly_mapping_output/mapped/{id}.bam",
    output:
        OUTDIR / "intermidiate_files/assembly_mapping_output/mapped_sorted/{id}.bam.sort",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", f"{str(OUTDIR)}/benchmark/") + "intermidiate_files_{id}_" + rulename
    log: 
        log=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_{id}_" + rulename,
        e=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_{id}_" + rulename+"_err",
        o=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_{id}_" + rulename+"_out"
    shell:
        """
    samtools sort --threads {threads} {input} -o {output} 2> {log.log}
    samtools index {output} 2>> {log.log}
    """

# extract coverage 
rule coverage:
    input:
        OUTDIR / "intermidiate_files/assembly_mapping_output/mapped_sorted/{id}.bam.sort",
    output:
        OUTDIR / "intermidiate_files/assembly_mapping_output/coverages/{id}_coverage.txt",
    threads: threads_fn("coverage")
    resources: walltime = walltime_fn("coverage"), mem_gb = mem_gb_fn("coverage")
    log:
        log=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_{id}_" + rulename,
        e=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_{id}_" + rulename+"_err",
        o=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_{id}_" + rulename+"_out"
    shell:
        "samtools coverage {input} > {output} 2> {log.log}}"

# 7. Run vamb to merge, split, and expand the hoods
rulename = "run_VAE"
rule run_VAE:
    input:
        contigs = "/home/projects/cu_10108/people/paupie/plasmaag_bangladesh_all_samples_wo_birth_1_2_2_2_based_on_louvain/results/candidate_plasmids_and_phage_contigs.fna",
        bamfiles = lambda wildcards: expand(OUTDIR / "intermidiate_files/assembly_mapping_output/mapped_sorted/{id}.bam.sort", id=sample_id["intermidiate_files"]),
        coverages = lambda wildcards: expand(OUTDIR / "intermidiate_files/assembly_mapping_output/coverages/{id}_coverage.txt", id=sample_id["intermidiate_files"])
    output:
        directory = directory(os.path.join(OUTDIR,"intermidiate_files", 'contrastive_VAE')),
        finished = os.path.join(OUTDIR,"intermidiate_files",'rule_completed_checks/run_VAE.finished'),
    params:
        cuda='--cuda' if CUDA else ''
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", f"{str(OUTDIR)}/benchmark/") + "intermidiate_files_" + rulename
    log: 
        log=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_" + rulename,
        e=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_" + rulename+"_err",
        o=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_" + rulename+"_out"
    shell:
        """
        rmdir {output.directory}
        {PLAMB_PRELOAD}
        vamb bin default --outdir {output.directory} --fasta {input.contigs} -p {threads} -e 2 -q 1 --bamfiles {input.bamfiles}\
        -m {MIN_CONTIG_LEN} {PLAMB_PARAMS}\
         {params.cuda} &> {log.log}
        touch {output.finished}
        """

