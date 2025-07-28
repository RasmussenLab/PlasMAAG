
# PlasMAAG 

PlasMAAG is a tool to recover **plasmids** and **organisms** from **metagenomic** samples, offering state-of-the-art plasmid reconstruction.
* On synthetic benchmark datasets, PlasMAAG reconstructs **50-121%** more near-complete confident plasmids than competing methods.
* On hospital sewage samples, PlasMAAG outperforms all other methods, reconstructing **33%** more plasmid sequences.
  
From FASTQ read files the tool will output:
* FASTA files of the most likely plasmids found in the reads.
* FASTA files of the most likely genomes found in the reads.  

_If you don't want PlasMAAG to assemble the reads, you can also pass in the paths to the asssemblies._

See the preprint for more information: ["Accurate plasmid reconstruction from metagenomics data using assembly-alignment graphs and contrastive learning"](https://www.biorxiv.org/content/10.1101/2025.02.26.640269v2.abstract)

## Quick Start :rocket:
Clone the repository and install the package using conda
```
git clone https://github.com/RasmussenLab/PlasMAAG
conda env create -n PlasMAAG --file=PlasMAAG/envs/PlasMAAG.yaml
```
To use the program activate the conda environment
```
conda activate PlasMAAG
```
To run the entire pipeline including assembly pass in a whitespace separated file containing the reads:
```
PlasMAAG --reads <read_file>  --output <output_directory> --threads <number_of_threads_to_use>
```
The <read_file> could look like:

``` 
read1                          read2
im/a/path/to/sample_1/read1    im/a/path/to/sample_1/read2
im/a/path/to/sample_2/read1    im/a/path/to/sample_2/read2
```
:heavy_exclamation_mark: Notice the header names are required to be: read1 and read2  

To dry run the pipeline before pass in the --dryrun flag

To run the pipeline from already assembled reads pass in a whitespace separated file containing the reads and the path to the spades assembly directories for each read pair.
```
PlasMAAG --reads_and_assembly_dir <reads_and_assembly_dir>  --output <output_directory> --threads <number_of_threads_to_use>
```
The `reads_and_assembly_dir` file could look like:
``` 
read1                          read2                         assembly_dir                                           
im/a/path/to/sample_1/read1    im/a/path/to/sample_1/read2   path/sample_1/Spades_output  
im/a/path/to/sample_2/read1    im/a/path/to/sample_2/read2   path/sample_2/Spades_output          
```
 :heavy_exclamation_mark: Notice the header names are required to be: read1, read2 and assembly_dir  

The path to the SPAdes output directory (under the assembly_dir column in the above example) must contain the following 3 files which Spades produces: 
| Description                         | File Name from Spades                               |
|:------------------------------------|:----------------------------------------|
| The assembled contigs               | `contigs.fasta`                         |
| The simplified assembly graphs      | `assembly_graph_after_simplification.gfa` |
| A metadata file                     | `contigs.paths`                         |

To see all options for the program run `PlasMAAG --help`

If interested in executing a testrun of PlasMAAG, please check the section ["Running the tool on test data
"](#Running-the-tool-on-test-data) or the Zenodo entry [here](https://zenodo.org/records/15263434). 


## Output files
The program produces three directories in the output directory choosen
```
< output directory >
├── intermidiate_files
├── log
└── results
```
The *results* directory contains:
````
results
├── candidate_plasmids.tsv : The candidate plasmids
├── candidate_genomes.tsv : The candidate chromosomes
└── scores.tsv : The aggregated scores for each cluster 
````
The `candidate_plasmids.tsv` and `candidate_genomes.tsv` files are formatted as:
````
clustername     contigname
nneighs_1051    Ssample0CNODE_198_length_19708_cov_59.381163
nneighs_1051    Ssample1CNODE_2317_length_2483_cov_58.855437
````
Here the sample names in the contignames (eg `sample0` in `Ssample0CNODE_198_length_19708_cov_59.381163`) refer to the order the reads were passed to the program. See example below:
``` 
read1                          read2
im/a/path/to/sample_1/read1    im/a/path/to/sample_1/read2  <-- Contigs from these reads would be called sample0
im/a/path/to/sample_2/read1    im/a/path/to/sample_2/read2  <-- Contigs from these reads would be called sample1
```

The *log* directory contains the output from the various rules called in the snakemake pipeline.
So for example `< output directory > /log/intermidiate_files_run_contrastive_VAE` would contain the log produced by running the `contrastive_VAE` rule in snakemake. 
 
The *intermidiate_files* directory contains the intermidiate files from the pipeline. 

## Running the tool on test data
Start of by setting up the tool:
```
git clone https://github.com/RasmussenLab/PlasMAAG
conda env create -n PlasMAAG --file=PlasMAAG/envs/PlasMAAG.yaml
conda activate PlasMAAG
```
Then download the test data
```
wget -O input_data.tar.gz .tar https://zenodo.org/records/15263434/files/input_data.tar.gz\?download\=1
tar -xvf input_data.tar.gz
```
Lastly run PlasMAAG. Here we pass in additional options to the VAE part of PlasMAAG (using <`vamb-arguments`), to have the tool work on the extremly small test dataset. Notice that the test dataset contains the `read_and_assembly_file.txt` file which is the configuration of the test data paths.
```
cd input_data
PlasMAAG --reads_and_assembly_dir read_and_assembly_file.txt --output test_run_PlasMAAG --threads 8 --vamb-arguments '-o C -e 200 -q 25 75 150 --seed 1'
```
Once the workflow finishes, several files and folders will be generated within the test_run_PlasMAAG directory. The final output files of the pipeline can be found in the the test_run_PlasMAAG/results directory, containing:
```
candidate_plasmids.tsv # The candidate plasmids
candidate_genomes.tsv # The candidate chromosomes
candidate_plasmids # Directory with the candidate plasmids fasta files
candidate_genomes # Directory with the candidate chromosomes fasta files
scores.tsv # The aggregated scores for each plasmid and genome cluster
```

## Advanced
### Using an already downloaded geNomad database
To use an already downloaded database, pass in a path to the genomad database with the ``` --genomad_db ``` argument

### Resources 

The pipeline can be configurated in: ``` config/config.yaml ```
Here, the resources for each rule can be configurated as follows
```
spades:
  walltime: "15-00:00:00"
  threads: 16
  mem_gb: 60
```
if no resourcess are configurated for a rule the defaults will be used which are also defined in: ``` config/config.yaml ```  as
```
default_walltime: "48:00:00"
default_threads: 16
default_mem_gb: 50
```
If these exceed the resourcess available they will be scaled down to match the hardware available. 

### Running on cluster [! needs to be tested also w.r.t to specific snakemake version as they changed how to submit to a cluster through snakemake !]
You can extend the arguments passed to snakemake by the '--snakemake_arguments' flag
This can then be used to have snakemake submit jobs to a cluster.
on PBS this could like:
```
PlasMAAG  <arguments> --snakemake_arguments \
    '--jobs 16 --max-jobs-per-second 5 --max-status-checks-per-second 5 --latency-wait 60 \
    --cluster "sbatch  --output={rule}.%j.o --error={rule}.%j.e \
    --time={resources.walltime} --job-name {rule} --cpus-per-task {threads} --mem {resources.mem_gb}G "'
```
### Running using snakemake CLI directly 
The pipeline can be run without using the CLI wrapper around snakemake. 
For using snakemake refer to the snakemake documentation: <https://snakemake.readthedocs.io/en/stable/>

#### Running from Reads using snakemake directly
To run the entire pipeline including assembly pass in a whitespace separated file containing the reads to snakemake using to the config flag in the snakemake CLI:
```
snakemake --use-conda --cores <number_of_cores> --snakefile <path_to_snakefile> --config read_file=<read_file>
```
The <read_file> could look like:

``` 
read1                          read2
im/a/path/to/sample_1/read1    im/a/path/to/sample_1/read2
im/a/path/to/sample_2/read1    im/a/path/to/sample_2/read2
```
:heavy_exclamation_mark: Notice the header names are required to be: read1 and read2  

#### Running from assembled reads using snakemake directly
To run the pipeline from allready assembled reads pass in a whitespace separated file containing the reads and the path to the spades assembly directories for each read pair to the config flag in snakemake.
```
snakemake --use-conda --cores <number_of_cores> --snakefile <path_to_snakefile> --config read_assembly_dir=<reads_and_assembly_dir_file>
```

The reads_and_assembly_dir_file could look like:
``` 
read1                          read2                         assembly_dir                                           
im/a/path/to/sample_1/read1    im/a/path/to/sample_1/read2   path/sample_1/Spades_output  
im/a/path/to/sample_2/read1    im/a/path/to/sample_2/read2   path/sample_2/Spades_output          
```
 :heavy_exclamation_mark: Notice the header names are required to be: read1, read2 and assembly_dir  

This assembly_dir directory filepath must contain the following 3 files which Spades produces: 
| Description                         | File Name from Spades                               |
|:------------------------------------|:----------------------------------------|
| The assembled contigs               | `contigs.fasta`                         |
| The simplified assembly graphs      | `assembly_graph_after_simplification.gfa` |
| A metadata file                     | `contigs.paths`                         |

#### Ressources for the different snakemake rules when using snakemake directly
To define resources for the specific snakemake rules either edit the snakemake.smk file directly or pass in a YAML config-file describing the resource use for each rule. 
```
snakemake <arguments> --configfile <config_file.yaml>
```
For a specification of the configfile refer to the ``` config/config.yaml ``` file.
Here the resources for each rule can be configurated as follows
```
spades:
  walltime: "15-00:00:00"
  threads: 16
  mem_gb: 60
```
if no resourcess are configurated for a rule the defaults will be used which are also defined in: ``` config/config.yaml ```  as
```
default_walltime: "48:00:00"
default_threads: 16
default_mem_gb: 50
```
#### Using an allready downloaded geNomad database 
To use an allready downloaded database, pass in a path to the genomad database using the config flag
```
snakemake <arguments> --config genomad_database=<path_to_genomad_database>
```

# TODO
- [ ] Lock in package versions for conda such that we do not need to solve the environments each time, and for reproducibility
- [ ] Test that cluster submit works
      - plasmid_bins/ *.fna 
      - nonplasmid_bins/ *.fna 
