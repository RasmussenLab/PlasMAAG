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
CONFIG_PATH = THIS_FILE_DIR / "config/config_shir.yaml"
# TODO should go to configfile
if CONFIG_PATH.exists():
    configfile: CONFIG_PATH


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


# Ensure consistent filename pattern for chunks
CHUNK_DIR = os.path.join(OUTDIR, "intermidiate_files/assembly_mapping_output/contig_chunks/")
CHUNK_PATTERN = os.path.join(CHUNK_DIR, "contigs_chunk_{id}.faa")


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
        os.path.join(OUTDIR, "intermidiate_files", "rule_completed_checks/fasta36/align_contigs.finished")


# chunk contigs.flt.fna.gz into chunks of size 100_000
rulename="chunk_contigs"
checkpoint chunk_contigs:
    input: "/home/projects/cu_10108/people/paupie/plasmaag_bangladesh_all_samples_wo_birth_1_2_2_2_based_on_louvain/results_new/prodigal/plasmids_filtered_proteins_renamed.faa"
    output: 
        directory(os.path.join(OUTDIR,"intermidiate_files/assembly_mapping_output/contig_chunks")),
        os.path.join(OUTDIR,"intermidiate_files",'rule_completed_checks/chunk_contigs.finished')
    params: os.path.join(OUTDIR, "intermidiate_files/assembly_mapping_output/contig_chunks/contigs_chunk_") 
    benchmark: config.get("benchmark", f"{str(OUTDIR)}/benchmark/") + "intermidiate_files" + rulename
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    log: 
        log=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_" + rulename,
        e=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_" + rulename+"_err",
        o=config.get("log", f"{str(OUTDIR)}/log/") + "intermidiate_files_" + rulename+"_out"
    shell:
        r"""
        mkdir -p {output[0]}
        awk -v prefix="{params}" -v chunk_size=50000 '
        /^>/ {{
            if (count % chunk_size == 0) {{
                if (out) close(out)
                out = sprintf("%s%04d.faa", prefix, ++file_id)
            }}
            count++
        }}
        {{ print | out }}
        ' {input} 2> {log.log}
        touch {output[1]}
        """

def get_chunk_ids(wildcards):
    ck = checkpoints.chunk_contigs.get(**wildcards)   # wait until chunking finished
    # Discover chunks by glob
    ids = glob_wildcards(CHUNK_PATTERN).id
    # Make deterministic order, so "upper triangle" is stable
    return sorted(ids)


def get_chunk_pairs_upper_with_diag(wc):
    ids = get_chunk_ids(wc)          # must return a sorted list of chunk ids, e.g. ['0001','0002',...]
    return [(a, b) for i, a in enumerate(ids) for b in ids[i:]]  # a <= b

def get_B_ids(wildcards):
    pairs = get_chunk_pairs(wildcards)
    return sorted({b for (a, b) in pairs})

rulename = "align_contigs_per_chunk"
rule align_contigs_per_chunk:
    input:
        # Query chunk for A
        lambda wildcards: os.path.join(
            CHUNK_DIR,
            f"contigs_chunk_{wildcards.a}.faa"
        ),
        lambda wildcards: os.path.join(
            CHUNK_DIR,
            f"contigs_chunk_{wildcards.b}.faa"
        ),
        os.path.join(OUTDIR,"intermidiate_files",'rule_completed_checks/chunk_contigs.finished')

    output:
        # BLAST output per chunk pair
        os.path.join(
            OUTDIR,
            "intermidiate_files",
            "fasta36",
            "chunk_pairwise",
            "fasta36_{a}_{b}.txt"
        ),
        # Finished flag
        os.path.join(
            OUTDIR,
            "intermidiate_files",
            "rule_completed_checks",
            "fasta36",
            "chunk_pairwise",
            "align_contigs_{a}_{b}.finished"
        )
    threads: threads_fn(rulename)
    resources:
        walltime = walltime_fn(rulename),
        mem_gb  = mem_gb_fn(rulename)
    benchmark:
        config.get("benchmark", f"{OUTDIR}/benchmark/") +
        "intermidiate_files_{a}_{b}_" + rulename
    log:
        log = config.get("log", f"{OUTDIR}/log/") +
              "intermidiate_files_{a}_{b}_" + rulename,
        e   = config.get("log", f"{OUTDIR}/log/") +
              "intermidiate_files_{a}_{b}_" + rulename + "_err",
        o   = config.get("log", f"{OUTDIR}/log/") +
              "intermidiate_files_{a}_{b}_" + rulename + "_out"
    shell:
        """
        module load tools 
        /home/projects/ku_00041/people/paupie/src/fasta36/bin/fasta36 -C 200  -m 8 {input[0]} {input[1]} | awk '$1 != $2' > {output[0]}  2>> {log.log}           
        touch {output[1]}
        """
rulename = "align_all_chunks"
rule align_all_chunks:
    input:
        # BLAST output files
        lambda wildcards: [
            os.path.join(
                OUTDIR,
                "intermidiate_files",
                "fasta36",
                "chunk_pairwise",
                f"fasta36_{a}_{b}.txt"
            )
            for a, b in get_chunk_pairs_upper_with_diag(wildcards)
        ],
        # Finished marker files
        lambda wildcards: [
            os.path.join(
                OUTDIR,
                "intermidiate_files",
                "rule_completed_checks/fasta36/chunk_pairwise",
                f"align_contigs_{a}_{b}.finished"
            )
            for a, b in get_chunk_pairs_upper_with_diag(wildcards)
        ]
    output:
        os.path.join(OUTDIR, "intermidiate_files", "fasta36", "fasta36_all_against_all.txt"),
        os.path.join(OUTDIR, "intermidiate_files", "rule_completed_checks/fasta36/align_contigs.finished")
    threads: threads_fn(rulename)
    resources:
        walltime = walltime_fn(rulename),
        mem_gb = mem_gb_fn(rulename)
    benchmark:
        config.get("benchmark", f"{OUTDIR}/benchmark/") + "intermidiate_files_" + rulename
    log:
        log = config.get("log", f"{OUTDIR}/log/") + "intermidiate_files_" + rulename,
        e   = config.get("log", f"{OUTDIR}/log/") + "intermidiate_files_" + rulename + "_err",
        o   = config.get("log", f"{OUTDIR}/log/") + "intermidiate_files_" + rulename + "_out"
    shell:
        """
        mkdir -p "$(dirname {output[0]})"
        : > "{output[0]}"

        find {OUTDIR}/intermidiate_files/fasta36/chunk_pairwise \
            -type f -name 'fasta36_*.txt' -print0 \
        | xargs -0 -r cat -- >> {output[0]} 2>> {log.log}

        touch {output[1]}
        """

# 3. Genereate nx graph from the alignment graph
rulename = "weighted_alignment_graph"
rule weighted_alignment_graph:
    input:
        os.path.join(OUTDIR,"intermidiate_files",'fasta36','fasta36_all_against_all.txt'),
        os.path.join(OUTDIR,"intermidiate_files",'rule_completed_checks/fasta36/align_contigs.finished')

