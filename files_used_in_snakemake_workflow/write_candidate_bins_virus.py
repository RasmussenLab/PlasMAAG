import argparse
import os

from Bio import SeqIO
from git_commit import get_git_commit
import gzip
import numpy as np 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write plasmid and non plasmid bins")

    # Add optional arguments with `--flags`
    parser.add_argument(
        "--cls_pl",
        required=True,
        type=str,
        help="Path to the plasmid clusters",
    )

    parser.add_argument(
        "--cls_org",
        required=True,
        type=str,
        help="Path to the organism clusters",
    )

    parser.add_argument(
        "--cls_vir",
        required=True,
        type=str,
        help="Path to the organism clusters",
    )

    parser.add_argument(
        "--contigs",
        required=True,
        type=str,
        help="Path to the fasta contigs file",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=str,
        help="Path to the results directory",
    )

    parser.add_argument(
        "--min_org_len",
        required=False,
        type=int,
        help="Min organism length",
        default=200_000,
    )

    parser.add_argument(
        "--min_vir_len",
        required=False,
        type=int,
        help="Min virus length",
        default=1,
    )

    parser.add_argument(
        "--min_plas_len",
        required=False,
        type=int,
        help="Min plasmid length",
        default=1,
    )

    
    ## Print git commit so we can debug
    commit_hash = get_git_commit(os.path.abspath(__file__))
    print("Git commit hash: " + commit_hash)

    # Parse the arguments
    args = parser.parse_args()
    print(args)

    # 
    plcl_cs_d = { cl:set() for cl,c in np.loadtxt(args.cls_pl,dtype=object,skiprows=1)}
    for cl,c in np.loadtxt(args.cls_pl,dtype=object,skiprows=1):
        plcl_cs_d[cl].add(c)
    plcl_len_d = { cl:np.sum([ int(c.split("length_")[1].split("_")[0]) for c in cs]) for cl,cs in plcl_cs_d.items() }
    c2plcl={ c:cl for cl,cs in plcl_cs_d.items() if plcl_len_d[cl] >= args.min_plas_len for c in cs}
    pl_cs = set(c2plcl.keys())
    
    
    orgcl_cs_d = { cl:set() for cl,c in np.loadtxt(args.cls_org,dtype=object,skiprows=1)}
    for cl,c in np.loadtxt(
        args.cls_org,dtype=object,skiprows=1):
        orgcl_cs_d[cl].add(c)
    
    orgcl_len_d = { cl:np.sum([ int(c.split("length_")[1].split("_")[0]) for c in cs]) for cl,cs in orgcl_cs_d.items() }
    c2orgcl={ c:cl for cl,cs in orgcl_cs_d.items() if orgcl_len_d[cl] >= args.min_org_len for c in cs}
    org_cs = set(c2orgcl.keys())

    vircl_cs_d = { cl:set() for cl,c in np.loadtxt(args.cls_vir,dtype=object,skiprows=1)}
    for cl,c in np.loadtxt(
        args.cls_vir,dtype=object,skiprows=1):
        vircl_cs_d[cl].add(c)
    
    vircl_len_d = { cl:np.sum([ int(c.split("length_")[1].split("_")[0]) for c in cs]) for cl,cs in vircl_cs_d.items() }
    c2vircl={ c:cl for cl,cs in vircl_cs_d.items() if vircl_len_d[cl] >= args.min_vir_len for c in cs}
    vir_cs = set(c2orgcl.keys())

    os.makedirs(os.path.join(args.outdir,"candidate_plasmids"))
    os.makedirs(os.path.join(args.outdir,"candidate_genomes"))
    os.makedirs(os.path.join(args.outdir,"candidate_viruses"))
    
    with gzip.open(args.contigs, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contig_name= record.id
            if contig_name not in org_cs.union(pl_cs).union(vir_cs):
                continue
            bin_dir =  "candidate_plasmids" if contig_name in pl_cs else "candidate_genomes" if contig_name in org_cs else "candidate_viruses"
            if bin_dir == "candidate_plasmids":
                bin_name= c2plcl[contig_name]
            elif bin_dir == "candidate_genomes":
                bin_name= c2orgcl[contig_name]
            else:
                bin_name= c2vircl[contig_name]
            bin_file = bin_name+".fna"
            bin_path = os.path.join(args.outdir,"%s/%s"%(bin_dir,bin_file))
            with open(bin_path, "a") as out:
                SeqIO.write(record, out, "fasta")

        
