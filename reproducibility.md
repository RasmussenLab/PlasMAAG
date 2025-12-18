Here are the instructions to be followed to reproduce the results presented on the PlasMAAG publication [https://www.biorxiv.org/content/10.1101/2025.02.26.640269v1]. Reproducibility steps are provided per figures.

# Figure 2.A,B,C,F.

Te generate the bins required for figure 2, we start with running PlasMAAG, and all the other binners, over the re-assembled CAMI2 datasets. 

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
3. Follow instructions to install PlasMAAG from zenodo [https://zenodo.org/records/17953597] 
4. Run PlasMAAG:
    ```bash
    conda activate PlasMAAG_zenodo # activate environment
    PlasMAAG --reads_and_assembly_dir read_and_assembly_dir_Airways.tsv  --output test_run_Airways -t 16 --vamb_arguments '-o C --seed 1 ’  # run plasmaag on the Airways dataset
    ```
    
5. PlasMAAG generates all the files necessary to run all the other binners:
    ```
    # contigs larger than 2kb across samples
   test_run_Airways/intermidiate_files/assembly_mapping_output/contigs.flt.fna.gz
    # sorted bam files
   test_run_Airways/intermidiate_files/assembly_mapping_output/mapped_sorted/*.bam.sort
    ```
7. Once we get the cluster files. We can evaluate the binning with BinBencher. Here [https://viralinstruction.com/BinBencherBackend.jl/dev/walkthrough/] is an extensive documentation about how to run BinBencher. In short, after installing julia and BinBencher, we just have to:
    ```bash
   # Download references
    wget https://www.erda.dk/archives/426ef81eb35ea07078bb0041ee186c84/references/Airways.json
    # Start and interactive julia session
    julia 
    ```
    ```julia
    using BinBencherBackend
    ref_f = 'Airways.json'
    ref_airways= open(i -> Reference(i), ref_f)
    Binning("test_runs_Airways/intermidiate_files/contrastive_VAE/vae_clusters_community_based_complete_and_circular_unsplit.tsv",ref_airways,binsplit_separator="C",filter_genomes=is_plasmid) 
    ```
        
   count the number of plasmid reconstructed from the community-based clustering, the number reported is the HQ bins, which should be somewhere around 310.

To evaluate the cellular genomes:
 ```julia
using BinBencherBackend
ref_f="Airways.json”
ref_airways= open(i -> Reference(i), ref_f)
    Binning("test_runs_Airways/intermidiate_files/contrastive_VAE/vae_clusters_density_unsplit.tsv",ref_airways,binsplit_separator="C",filter_genomes=is_organism)
```

count the number of cellular genomes reconstructed from the density-based clustering, the number reported is the HQ bins, which should be somewhere around 53.

        
10. To run SCAPP, first you have to install it [https://github.com/Shamir-Lab/SCAPP?tab=readme-ov-file#installation], and then you are ready to run it:
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
        
12. Evaluate SCAPP binning with BinBencher
    ```bash
    # Download references
    wget https://www.erda.dk/archives/426ef81eb35ea07078bb0041ee186c84/references/scapp_references/Airways.json
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

Te generate the bins required for figure 4.B, we follow a similar procedure than with Figure 2. We start with running PlasMAAG, and all the other binners, over the re-assembled CAMI2 datasets. 

1. Download input files, reads, and assembly directories.
    1. ```bash
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
   ```
   # contigs larger than 2kb across samples
   test_run_DARWIN/intermidiate_files/assembly_mapping_output/contigs.flt.fna.gz
   # sorted bam files
   test_run_DARWIN/intermidiate_files/assembly_mapping_output/mapped_sorted/*.bam.sort 
   ```
6. Once we get the cluster files. We can evaluate the binning with BinBencher. Here [https://viralinstruction.com/BinBencherBackend.jl/dev/walkthrough/] is an extensive documentation about how to run BinBencher. In short, after installing julia and BinBencher, we just have to:
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
        
7. To run SCAPP, first you have to install it [https://github.com/Shamir-Lab/SCAPP?tab=readme-ov-file#installation], and then you are ready to run it:        
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
        
9. Evaluate SCAPP binning with BinBencher
    ```bash
   # Download references
    wget https://www.erda.dk/archives/426ef81eb35ea07078bb0041ee186c84/references/scapp_references/DARWIN_long_reads_scapp.json
    # Start and interactive julia session
    julia 
    ```                
    ```julia
    using BinBencherBackend
    ref_f="DARWIN_long_reads_scapp.json"
    ref_darwin_scapp= open(i -> Reference(i), ref_f)
    Binning("scapp_DARWIN/cycles_confident_clusters.tsv",ref_darwin_scapp,binsplit_separator="C",filter_genomes=is_plasmid) 
    ```
   count the number of plasmid reconstructed from the SCAPP confident, the number reported is the HQ bins, which should be somewhere around 17.
