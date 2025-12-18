Here are the instructions to be followed to reproduce the results presented on the PlasMAAG publication [https://www.biorxiv.org/content/10.1101/2025.02.26.640269v1]. Reproducibility steps are provided per figures.

# Figure 2.A,B,C,F.

To generate the bins required for figure 2, we start with running PlasMAAG, and all the other binners, over the re-assembled CAMI2 datasets.

1. Download the input files, reads, and assembly directories.

```bash
# download assembly directories and decompress them
mkdir assemblies; cd assemblies
wget https://www.erda.dk/archives/73d8780489ac5364d8cc4e8093bc9303/CAMI2_assemblies_and_PlasMAAG_input/Airways.tar.gz
tar -xvf Airways.tar.gz; cd ..
# download reads directories and decompress them
mkdir reads; cd reads
wget https://erda.ku.dk/archives/826fe4d8889f88db2ec20058f9eaa015/reassembled_CAMI_reads.tar.gz
tar -xvf reassembled_CAMI_reads.tar.gz
cd ..
# download PlasMAAG input
wget https://www.erda.dk/archives/73d8780489ac5364d8cc4e8093bc9303/CAMI2_assemblies_and_PlasMAAG_input/read_and_assembly_dir_Airways.tsv
```

2. Follow instructions to install PlasMAAG from zenodo [https://zenodo.org/records/17953597]

3. Run PlasMAAG:

```bash
conda activate PlasMAAG_zenodo # activate environment
PlasMAAG --reads_and_assembly_dir read_and_assembly_dir_Airways.tsv  --output test_run_Airways -t 16 --vamb_arguments '-o C --seed 1 ’  # run plasmaag on the Airways dataset
```

4. PlasMAAG generates all the files necessary to run all the other binners:

```bash
# contigs larger than 2kb across samples
test_run_Airways/intermidiate_files/assembly_mapping_output/contigs.flt.fna.gz
# sorted bam files
test_run_Airways/intermidiate_files/assembly_mapping_output/mapped_sorted/*.bam.sort
```

5. Once we get the cluster files. We can evaluate the binning with BinBencher.
   The BinBencher documentation can be found here [https://viralinstruction.com/BinBencherBackend.jl/v0.3.4/], and the BinBencher installation instructions can be found here [https://github.com/jakobnissen/BinBencher.jl].
   Once it's installed, we can run it from command line, like so:

```bash
# Download references
wget https://www.erda.dk/archives/426ef81eb35ea07078bb0041ee186c84/references/Airways.json
# run BinBencher
binbench bench -s C --keep-flags plasmid out1 Airways.json test_run_Airways/intermidiate_files/contrastive_VAE/vae_clusters_community_based_complete_and_circular_unsplit.tsv
```

The result is the fifth element (precision: 0.95) of the fourth element (recall: 0.9) of the first element (genome level) of the "genomes_genomic_recall" field in the file `recovery.json`.
One way to extract this number from command line is with the JSON reading tool `jq`.
The number reported is the HQ plasmids reconstructed from community based clustering, which should be somewhere around 310.
```bash
jq '.genomes_genomic_recall[0][3][4]' out1/recovery.json
```
We can also do it in an interactive session:
 ```bash
 # Start and interactive julia session
 julia 
 ```
 ```julia
 using BinBencherBackend
 ref_f = "Airways.json"
 ref_airways= open(i -> Reference(i), ref_f)
 Binning("test_run_Airways/intermidiate_files/contrastive_VAE/vae_clusters_community_based_complete_and_circular_unsplit.tsv",ref_airways,binsplit_separator="C",filter_genomes=is_plasmid) 
 ```

To evaluate the cellular genomes, as found using density based clustering:

```bash
binbench bench -s C --keep-flags organism out2 Airways.json test_run_Airways/intermidiate_files/contrastive_VAE/vae_clusters_density_unsplit.tsv
jq '.genomes_genomic_recall[0][3][4]' out2/recovery.json
```
We can also do it in an interactive session:
 ```bash
 # Start and interactive julia session
 julia 
 ```
 ```julia
using BinBencherBackend
ref_f = "Airways.json"
ref_airways= open(i -> Reference(i), ref_f)
Binning("test_run_Airways/intermidiate_files/contrastive_VAE/vae_clusters_density_unsplit.tsv",ref_airways,binsplit_separator="C",filter_genomes=is_organism)
```
The answer should be around 53.

6. To run SCAPP, first you have to install it [https://github.com/Shamir-Lab/SCAPP?tab=readme-ov-file#installation], and then you are ready to run it:

```bash
# Run SCAPP for the re-assembled CAMI2 Airways dataset
dir_asm=assemblies/Airways/
mkdir scapp_Airways # create the directory for the Airways samples
for s in $(ls $dir_asm)
do
graph_file=${dir_asm}/"$s"/assembly_graph.fastg
outdir=scapp_Airways/"$s"
only_number_s=$(echo $s | sed 's=S==g')
reads1=reads/Airways/simulation/*_sample_${only_number_s}/reads/reads_noninterlaced/anonymous_reads_clean_1.fq
reads2=reads/Airways/simulation/*_sample_${only_number_s}/reads/reads_noninterlaced/anonymous_reads_clean_2.fq

scapp -g $graph_file -o $outdir  -r1 $reads1 -r2 $reads2 -p 16

## aggregate scapp results across samples into a cluster.tsv file, which we can then benchmark with BinBencher
clusters_cycles=scapp_Airways/cycles_clusters.tsv
clusters_confident_cycles=scapp_Airways/cycles_confident_clusters.tsv
echo -e "clustername\tcontigname" > $clusters_cycles
for s in $(ls scapp_Airways/S*)
do
# cycles
cycles_fasta=scapp_Airways/"$s"/intermediate_files/assembly_graph.cycs.fasta
grep '^>' $cycles_fasta *|* sed "s=>="$s"=g" | awk '{print $1 "\t" $1}' >> $clusters_cycles
# confident cycles
confident_cycles_fasta=scapp_Airways/"$s"/assembly_graph.confident_cycs.fasta
grep '^>' $confident_cycles_fasta *|* sed "s=>="$s"=g" | awk '{print $1 "\t" $1}' >> $clusters_confident_cycles
done
```
7. Evaluate SCAPP binning with BinBencher
```bash
# Download references
wget https://www.erda.dk/archives/426ef81eb35ea07078bb0041ee186c84/references/scapp_references/Airways.json

# Check the number of plasmid reconstructed from the SCAPP confident, which should be somewhere around 189.
binbench bench -s C --keep-flags plasmid out3 Airways.json scapp_Airways/cycles_confident_clusters.tsv
jq '.genomes_genomic_recall[0][3][4]' out3/recovery.json

# To evaluate SCAPP cycles, which should be around 254:
binbench bench -s C --keep-flags plasmid out4 Airways.json scapp_Airways/cycles_clusters.ts
jq '.genomes_genomic_recall[0][3][4]' out4/recovery.json
```
Or in an interactive session:
 ```bash
 # Start and interactive julia session
 julia 
 ```        
 ```julia
 using BinBencherBackend
 ref_f="Airways.json”
 ref_airways_scapp= open(i -> Reference(i), ref_f)
 Binning("scapp_Airways/cycles_confident_clusters.tsv",ref_airways_scapp,binsplit_separator="C",filter_genomes=is_plasmid) 
 ```
 count the number of plasmid reconstructed from the SCAPP confident, the number reported is the HQ bins, which should be somewhere around 189.

 To evaluate SCAPP_cycles:
 ```julia
 using BinBencherBackend
 ref_f =  'Airways.json'
 ref_airways_scapp= open(i -> Reference(i), ref_f)
 Binning("scapp_Airways/cycles_clusters.tsv",ref_airways_scapp,binsplit_separator="C",filter_genomes=is_plasmid) 
 ```
 count the number of plasmid reconstructed from the SCAPP cycles, the number reported is the HQ bins, which should be somewhere around 254.

# Figure 4.B

To generate the bins required for figure 4.B, we follow a similar procedure than with Figure 2. We start with running PlasMAAG, and all the other binners, over the re-assembled CAMI2 datasets.

1. Download input files, reads, and assembly directories.
```bash
# download assembly, reads, and PlasMAAG input file
wget https://www.erda.dk/archives/753b8c039aa18adc2956973d376de97f/DARWIN_assemblies_reads_and_PlasMAAG_input/assemblies.tar.gz
tar -xvf assemblies.tar.gz
wget https://www.erda.dk/archives/753b8c039aa18adc2956973d376de97f/DARWIN_assemblies_reads_and_PlasMAAG_input/reads.tar.gz
tar -xvf reads.tar.gz
wget https://www.erda.dk/archives/753b8c039aa18adc2956973d376de97f/DARWIN_assemblies_reads_and_PlasMAAG_input/read_and_assembly_dir.tsv
```

2. Follow instructions to install PlasMAAG from zenodo [https://zenodo.org/records/17953597]

3. Run PlasMAAG:

```bash
conda activate PlasMAAG_zenodo # activate environment
PlasMAAG --reads_and_assembly_dir read_and_assembly_dir.tsv  --output test_run_DARWIN -t 16 --vamb_arguments '-o C --seed 1 ’
```

4. PlasMAAG generates all the files necessary to run all the other binners:

```bash
# contigs larger than 2kb across samples
test_run_DARWIN/intermidiate_files/assembly_mapping_output/contigs.flt.fna.gz
# sorted bam files
test_run_DARWIN/intermidiate_files/assembly_mapping_output/mapped_sorted/*.bam.sort
```

5. Once we get the cluster files. We can evaluate the binning with BinBencher.
   See the instructions above for how to install BinBencher, and its documentation.
   After installing, run:

```bash
# Download references
wget https://www.erda.dk/archives/426ef81eb35ea07078bb0041ee186c84/references/DARWIN_long_reads.json

# Benchmark number of long-read assemblies reconstructed from the community-based clustering
binbench bench -s C out5 DARWIN_long_reads.json test_runs_DARWIN/intermidiate_files/contrastive_VAE/vae_clusters_community_based_complete_and_circular_unsplit.tsv
# get the number from the json file with jq
# This number should be around 63
jq '.genomes_genomic_recall[0][3][4]' out5/recovery.json
```

6. To run SCAPP, first you have to install it [https://github.com/Shamir-Lab/SCAPP?tab=readme-ov-file#installation], and then you are ready to run it:

```bash
# Run SCAPP for the DARWIN dataset
dir_asm=assemblies
mkdir scapp_DARWIN # create the directory for the DARWIN samples
for s in $(ls $dir_asm)
do
graph_file=${dir_asm}/"$s"/assembly_graph.fastg
outdir=scapp_DARWIN/"$s"
sample_name_without_prefix=$(echo $s | cut -f 2 -d C)
reads1=reads/"$sample_name_without_prefix"_1.qc.fastq.gz
reads2=reads/"$sample_name_without_prefix"_2.qc.fastq.gz

# run scapp
scapp -g $graph_file -o $outdir  -r1 $reads1 -r2 $reads2 -p 8

## aggregate scapp results across samples into a cluster.tsv file, which we can then benchmark with BinBencher
clusters_cycles=scapp_DARWIN/cycles_clusters.tsv
clusters_confident_cycles=scapp_DARWIN/cycles_confident_clusters.tsv
echo -e "clustername\tcontigname" > $clusters_cycles
for s in $(ls scapp_DARWIN/S*)
do
# cycles
cycles_fasta=scapp_DARWIN/"$s"/intermediate_files/assembly_graph.cycs.fasta
grep '^>' $cycles_fasta *|* sed "s=>="$s"=g" | awk '{print $1 "\t" $1}' >> $clusters_cycles
# confident cycles
confident_cycles_fasta=scapp_DARWIN/"$s"/assembly_graph.confident_cycs.fasta
grep '^>' $confident_cycles_fasta *|* sed "s=>="$s"=g" | awk '{print $1 "\t" $1}' >> $clusters_confident_cycles
done
```

7. Evaluate SCAPP binning with BinBencher
```bash
# Download references
wget https://www.erda.dk/archives/426ef81eb35ea07078bb0041ee186c84/references/scapp_references/DARWIN_long_reads_scapp.json


# Benchmark number of plasmid reconstructed from the SCAPP confident
binbench bench -s C out6 DARWIN_long_reads_scapp.json scapp_DARWIN/cycles_confident_clusters.tsv
# get the number from the json file with jq
# This number should be around 17
jq '.genomes_genomic_recall[0][3][4]' out6/recovery.json
```
Or in an interactive session:
 ```bash
# Download references
 wget https://www.erda.dk/archives/426ef81eb35ea07078bb0041ee186c84/references/DARWIN_long_reads.json
 # Start and interactive julia session
 julia 
 ```        
 ```julia
 using BinBencherBackend
 ref_f="DARWIN_long_reads.json"
 ref_darwin= open(i -> Reference(i), ref_f)
Binning("test_runs_DARWIN/intermidiate_files/contrastive_VAE/vae_clusters_community_based_complete_and_circular_unsplit.tsv",ref_darwin,binsplit_separator="C") 
 ```
count the number of long-read assemblies reconstructed from the community-based clustering, the number reported is the HQ bins, which should be somewhere around 63.

# Generating your own BinBencher reference files
The examples above show how to make use of pre-defined BinBencher reference JSON files.
To run and benchmark PlasMAAG on your own samples, you may want to create your own references.

!!! note
    BinBencher reference files store a _ground truth reference_ for a dataset,
    and can therefore only be used if you have some sort of ground truth, e.g.
    a simulated data, or in the case of PlasMAAG, paired long/short read. 

You will need:
* The set of genomes in your sample, as FASTA files
* The taxonomy of these genomes, in a tree, if you want benchmarking on other phylogenetic levels than genome
* The set of sequences used for binning as a FASTA file
* The mapping positions of sequence to the genomes

For specifics about how to create your own reference JSON, see the BinBencher documentation [https://viralinstruction.com/BinBencherBackend.jl/v0.3.4/make_ref/]
