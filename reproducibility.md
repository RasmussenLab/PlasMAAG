Here are the instructions to be followed to reproduce the results presented on the PlasMAAG publication [https://www.biorxiv.org/content/10.1101/2025.02.26.640269v1]. Reproducibility steps are provided per figures.

# Figure 2.A,B,C,F.

Te generate the bins required for figure 2, we start with running PlasMAAG, and all the other binners, over the re-assembled CAMI2 datasets. 

1. Download the reads, and assembly directories.
    1. `wget [ERDA link to CAMI2 assembly dirs]; tar -xvf spades_output.tar.gz # download assembly directories and decompress them`
    2. `wget [ERDA link to CAMI2 reads]; tar -xvf # download reads directories and decompress them `
2. Download PlasMAAG input file:
    1. `wget [ERDA linkt to CAMI2 plasmaag input files]`
3. Follow instructions to install PlasMAAG from zenodo [https://zenodo.org/records/17953597] 
4. Run PlasMAAG:
    1. `conda activate PlasMAAG_zenodo # activate environment `
    2. `PlasMAAG --reads_and_assembly_dir read_and_assembly_dir_Airways.tsv  --output test_run_Airways -t 16 --vamb_arguments '-o C --seed 1’ --genomad_thr '0.1’  # run plasmaag on the Airways dataset `
5. PlasMAAG generates all the files necessary to run all the other binners:
    1. `test_run_Airways/intermidiate_files/assembly_mapping_output/contigs.flt.fna.gz # contigs larger than 2kb across samples `
    2. `test_run_Airways/intermidiate_files/assembly_mapping_output/mapped_sorted/*.bam.sort # sorted bam files `
6. Once we get the cluster files. We can evaluate the binning with BinBencher. Here [https://viralinstruction.com/BinBencherBackend.jl/dev/walkthrough/] is an extensive documentation about how to run BinBencher. In short, after installing julia and BinBencher, we just have to:
    1. `wget [ERDA CAMI2 references] # Download references`
    2. `julia # Start and interactive julia session`
        
        ```julia
        using BinBencherBackend
        ref_f =  “references/Airways.json”
        ref_airways= open(i -> Reference(i), ref_f)
        Binning("test_runs/Airways/intermidiate_files/contrastive_VAE/vae_clusters_community_based_complete_and_circular_unsplit.tsv",refa,binsplit_separator="C",filter_genomes=is_plasmid) 
        ```
        
        count the number of plasmid reconstructed from the community-based clustering, the number reported is the HQ bins, which should be somewhere around 310.
        
7. To run SCAPP, first you have to install it [https://github.com/Shamir-Lab/SCAPP?tab=readme-ov-file#installation], and then you are ready to run it:
    1. Run SCAPP for the re-assembled CAMI2 dataset
        
        ```bash
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
        
8. Evaluate SCAPP binning with BinBencher
    1. `wget [ERDA CAMI2 SCAPP references] # Download references`
    2. `julia # Start and interactive julia session`
        
        ```julia
        using BinBencherBackend
        ref_f =  “references_scapp/Airways.json”
        ref_airways= open(i -> Reference(i), ref_f)
        Binning("scapp_Airways/cycles_confident_clusters.tsv",refa,binsplit_separator="C",filter_genomes=is_plasmid) 
        ```
        
        count the number of plasmid reconstructed from the SCAPP confident, the number reported is the HQ bins, which should be somewhere around 189.
        

# Figure 3.B

Te generate the bins required for figure 3.B, we follow a similar procedure than with Figure 2. We start with running PlasMAAG, and all the other binners, over the re-assembled CAMI2 datasets. 

1. Download the reads, and assembly directories.
    1. `wget [ERDA link to DARWIN assembly dirs]; tar -xvf spades_output.tar.gz    # download assembly directories and decompress them`
    2. `wget [ERDA link to DARWIN reads]; tar -xvf # download reads directories and decompress them`
2. Download PlasMAAG input file:
    1. `wget [ERDA linkt to DARWIN plasmaag input files]`
3. Follow instructions to install PlasMAAG from zenodo [https://zenodo.org/records/17953597] 
4. Run PlasMAAG:
    1. `conda activate PlasMAAG_zenodo # activate environment`
    2. `PlasMAAG --reads_and_assembly_dir read_and_assembly_dir_DARWIN.tsv  --output test_run_DARWIN -t 16 --vamb_arguments '-o C --seed 1’ --genomad_thr '0.95’  `
5. PlasMAAG generates all the files necessary to run all the other binners:
    1. `test_run_DARWIN/intermidiate_files/assembly_mapping_output/contigs.flt.fna.gz # contigs larger than 2kb across samples`
    2. `test_run_DARWIN/intermidiate_files/assembly_mapping_output/mapped_sorted/*.bam.sort # sorted bam files`
6. Once we get the cluster files. We can evaluate the binning with BinBencher. Here [https://viralinstruction.com/BinBencherBackend.jl/dev/walkthrough/] is an extensive documentation about how to run BinBencher. In short, after installing julia and BinBencher, we just have to:
    1. `wget [ERDA DARWIN references] # Download references`
    2. `julia # Start and interactive julia session`
        
        ```julia
        using BinBencherBackend
        ref_f =  “references/DARWIN_long_reads.json”
        ref_darwin= open(i -> Reference(i), ref_f)
        Binning("test_runs_DARWIN/intermidiate_files/contrastive_VAE/vae_clusters_community_based_complete_and_circular_unsplit.tsv",ref_darwin,binsplit_separator="C") 
        ```
        
        count the number of long-read assemblies reconstructed from the community-based clustering, the number reported is the HQ bins, which should be somewhere around 62.
        
7. To run SCAPP, first you have to install it [https://github.com/Shamir-Lab/SCAPP?tab=readme-ov-file#installation], and then you are ready to run it:
    1. Run SCAPP for the re-assembled DARWIN dataset
        
        ```bash
        dir_asm=assemblies/DARWIN/
        mkdir scapp_DARWIN # create the directory for the DARWIN samples
        for s in $(ls $dir_asm)
        do
        graph_file=${dir_asm}/"$s"/assembly_graph.fastg
        outdir=scapp_DARWIN/"$s"
        only_number_s=$(echo $s | sed 's=S==g')
        reads1=reads/DARWIN/simulation/*_sample_${only_number_s}/reads/reads_noninterlaced/anonymous_reads_clean_1.fq
        reads2=reads/DARWIN/simulation/*_sample_${only_number_s}/reads/reads_noninterlaced/anonymous_reads_clean_2.fq
        
        scapp -g $graph_file -o $outdir  -r1 $reads1 -r2 $reads2 -p 16
        
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
        
8. Evaluate SCAPP binning with BinBencher
    1. `wget [ERDA CAMI2 SCAPP references] # Download references`
    2. `julia # Start and interactive julia session`
        
        ```julia
        using BinBencherBackend
        ref_f =  “references_scapp/DARWIN_long_reads.json”
        ref_DARWIN= open(i -> Reference(i), ref_f)
        Binning("scapp_DARWIN/cycles_confident_clusters.tsv",ref_DARWIN,binsplit_separator="C",filter_genomes=is_plasmid) 
        ```
        
        count the number of plasmid reconstructed from the SCAPP confident, the number reported is the HQ bins, which should be somewhere around 17.
