import argparse
import os

from Bio import SeqIO
from git_commit import get_git_commit
import gzip
import numpy as np 
from pathlib import Path
from collections import OrderedDict



def split_clusters_by_sample(clusters):
    cls_split = dict()
    for i, (cl_id, cl_cs) in enumerate(clusters.items()):
        samples_in_cl = set([c.split("C")[0] for c in cl_cs])

        cl_S_d = {S: set() for S in samples_in_cl}

        for c in cl_cs:
            cl_S_d[c.split("C")[0]].add(c)

        for S, cs in cl_S_d.items():
            cls_split["%s_%s" % (str(cl_id), S)] = cs

    return cls_split

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
    plcl_cs_d_ = { cl:set() for cl,c in np.loadtxt(args.cls_pl,dtype=object,skiprows=1)}
    for cl,c in np.loadtxt(args.cls_pl,dtype=object,skiprows=1):
        plcl_cs_d_[cl].add(c)
    plcl_cs_d = split_clusters_by_sample(plcl_cs_d_)
    plcl_len_d = { cl:np.sum([ int(c.split("length_")[1].split("_")[0]) for c in cs]) for cl,cs in plcl_cs_d.items() }
    c2plcl={ c:cl for cl,cs in plcl_cs_d.items() if plcl_len_d[cl] >= args.min_plas_len for c in cs}
    pl_cs = set(c2plcl.keys())
    
    
    orgcl_cs_d_ = { cl:set() for cl,c in np.loadtxt(args.cls_org,dtype=object,skiprows=1)}
    for cl,c in np.loadtxt(
        args.cls_org,dtype=object,skiprows=1):
        orgcl_cs_d_[cl].add(c)
    orgcl_cs_d = split_clusters_by_sample(orgcl_cs_d_)
    orgcl_len_d = { cl:np.sum([ int(c.split("length_")[1].split("_")[0]) for c in cs]) for cl,cs in orgcl_cs_d.items() }
    c2orgcl={ c:cl for cl,cs in orgcl_cs_d.items() if orgcl_len_d[cl] >= args.min_org_len for c in cs}
    org_cs = set(c2orgcl.keys())

    vircl_cs_d_ = { cl:set() for cl,c in np.loadtxt(args.cls_vir,dtype=object,skiprows=1)}
    for cl,c in np.loadtxt(
        args.cls_vir,dtype=object,skiprows=1):
        vircl_cs_d_[cl].add(c)
    vircl_cs_d = split_clusters_by_sample(vircl_cs_d_)
    vircl_len_d = { cl:np.sum([ int(c.split("length_")[1].split("_")[0]) for c in cs]) for cl,cs in vircl_cs_d.items() }
    c2vircl={ c:cl for cl,cs in vircl_cs_d.items() if vircl_len_d[cl] >= args.min_vir_len for c in cs}
    vir_cs = set(c2vircl.keys())


    # Precompute a single lookup: contig_name -> (bin_dir, bin_name)
    contig_to_bin = {}

    # Fill in plasmids
    for contig_name, bin_name in c2plcl.items():
        contig_to_bin[contig_name] = ("candidate_plasmids", bin_name)

    # Fill in genomes
    for contig_name, bin_name in c2orgcl.items():
        contig_to_bin[contig_name] = ("candidate_genomes", bin_name)

    # Fill in viruses
    for contig_name, bin_name in c2vircl.items():
        contig_to_bin[contig_name] = ("candidate_viruses", bin_name)

    outdir = Path(args.outdir)

    # Ensure directories exist (no error if already present)
    for d in ("candidate_plasmids", "candidate_genomes", "candidate_viruses"):
        (outdir / d).mkdir(parents=True, exist_ok=True)



    # LRU cache for file handles
    MAX_OPEN_FILES = 512  # tune based on your environment; keep well below OS limit
    handles = OrderedDict()  # key: Path, value: open file object

    def get_handle(path: Path):
        """Get a file handle with LRU eviction when exceeding MAX_OPEN_FILES."""
        # If handle exists, mark it as recently used
        fh = handles.get(path)
        if fh is not None:
            # move to the end to denote recent use
            handles.move_to_end(path, last=True)
            return fh

        # Need to open a new handle
        fh = path.open("a")
        handles[path] = fh
        handles.move_to_end(path, last=True)

        # Evict least-recently-used if over capacity
        if len(handles) > MAX_OPEN_FILES:
            old_path, old_fh = handles.popitem(last=False)
            try:
                old_fh.flush()
            finally:
                old_fh.close()
        return fh

    try:
        with gzip.open(args.contigs, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                info = contig_to_bin.get(record.id)
                if info is None:
                    continue
                bin_dir, bin_name = info
                bin_path = outdir / bin_dir / f"{bin_name}.fna"
                out_fh = get_handle(bin_path)
                SeqIO.write(record, out_fh, "fasta")
    finally:
        # Close any remaining open handles
        for fh in handles.values():
            try:
                fh.flush()
            finally:
                fh.close()


    # os.makedirs(os.path.join(args.outdir,"candidate_plasmids"))
    # os.makedirs(os.path.join(args.outdir,"candidate_genomes"))
    # os.makedirs(os.path.join(args.outdir,"candidate_viruses"))

    # with gzip.open(args.contigs, "rt") as handle:
    #     for record in SeqIO.parse(handle, "fasta"):
    #         contig_name= record.id
    #         if contig_name not in org_cs.union(pl_cs).union(vir_cs):
    #             continue
    #         bin_dir =  "candidate_plasmids" if contig_name in pl_cs else "candidate_genomes" if contig_name in org_cs else "candidate_viruses"
    #         if bin_dir == "candidate_plasmids":
    #             bin_name= c2plcl[contig_name]
    #         elif bin_dir == "candidate_genomes":
    #             bin_name= c2orgcl[contig_name]
    #         else:
    #             bin_name= c2vircl[contig_name]
    #         bin_file = bin_name+".fna"
    #         bin_path = os.path.join(args.outdir,"%s/%s"%(bin_dir,bin_file))
    #         with open(bin_path, "a") as out:
    #             SeqIO.write(record, out, "fasta")

        
