
# PlasMAAG 

PlasMAAG is a tool to recover **plasmids** and **organisms** from **metagenomic** samples, offering state-of-the-art plasmid reconstruction.
* On synthetic benchmark datasets, PlasMAAG reconstructs **50-121%** more near-complete confident plasmids than competing methods.
* On hospital sewage samples, PlasMAAG outperforms all other methods, reconstructing **33%** more plasmid sequences.
  
From FASTQ read files the tool will output:
* FASTA files of the most likely plasmids found in the reads (For genomes these are refered to as MAG's).
* FASTA files of the most likely genomes found in the reads (Normally refered to as MAG's).  

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
├── candidate_plasmids/ : A directory containing the the candidate plasmids
├── candidate_genomes/ : A directory containing the candidate chromosomes
├── candidate_plasmids.tsv : An overview of which contigs are binned together as candidate plasmids
├── candidate_genomes.tsv : An overview of which contigs are binned together as candidate chromosomes
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
wget -O input_data.tar.gz https://zenodo.org/records/15263434/files/input_data.tar.gz\?download\=1
tar -xvf input_data.tar.gz
```
Lastly run PlasMAAG.  
_Here we pass in additional options to the VAE part of PlasMAAG (using <`--vamb_arguments`), to have the tool work on the extremly small test dataset._
```
cd input_data
PlasMAAG --reads_and_assembly_dir read_and_assembly_file.txt --output test_run_PlasMAAG --threads 8 --vamb_arguments '-o C -e 200 -q 25 75 150 --seed 1'
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

### Running using snakemake CLI directly 
The pipeline can be run without using the CLI wrapper around snakemake. 
For using snakemake refer to the snakemake documentation: <https://snakemake.readthedocs.io/en/stable/>

#### Running from Reads using snakemake directly
To run the entire pipeline including assembly pass in a whitespace separated file containing the reads to snakemake using to the config flag in the snakemake CLI:
```
snakemake --use-conda --cores <number_of_cores> --snakefile <path_to_snakefile> --config read_file=<read_file> output_directory=<output_directory>
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
snakemake --use-conda --cores <number_of_cores> --snakefile <path_to_snakefile> --config read_assembly_dir=<reads_and_assembly_dir_file>  output_directory=<output_directory>
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

### Running on a cluster with snakemake submiting jobs 
For running PlasMAAG on a cluster with snakemake submiting jobs see the documentation for snakemake [here](https://snakemake.readthedocs.io/en/v7.19.1/executing/cluster.html)  
An example is provided below for reference using slurm running PlasMAAG from reads:
Start off by installing the cluster-generic executor plugin for snakemake
```
pip install snakemake-executor-plugin-cluster-generic
```
Then run the PlasMAAG snakemake pipeline:
```
snakemake --use-conda --snakefile <path_to_snakefile> --config read_assembly_dir=<reads_and_assembly_dir_file> output_directory=<output_directory> \
  --jobs 2 --max-jobs-per-second 5 --max-status-checks-per-second 5 --latency-wait 60 \
  --executor cluster-generic --cluster-generic-submit-cmd 'sbatch --job-name {rule} --time={resources.walltime} --cpus-per-task {threads} --mem {resources.mem_gb}G'
```
#### Ressources for the different snakemake rules when using snakemake directly
To define resources for the specific snakemake rules edit the `config/config.yaml` file
For more information see the ["Resources" section](#Resources).

#### Using an allready downloaded geNomad database 
To use an allready downloaded database, pass in a path to the genomad database using the config flag
```
snakemake <arguments> --config genomad_database=<path_to_genomad_database>
```

