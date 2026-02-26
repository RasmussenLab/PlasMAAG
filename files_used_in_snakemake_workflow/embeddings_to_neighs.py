import argparse
import os
import pickle

# from scipy.spatial.distance import cdist
import time

import networkx as nx
import numpy as np
import torch
from git_commit import get_git_commit

import vamb

from collections import defaultdict

# def find_neighbours_optimized(
#     embeddings_bincontigs,
#     contignames,
#     ccs_graph_d,
#     radius_clustering=0.2,
#     build_graph=False,
# ):
#     """Optimized neighbor finder that only compares contigs belonging 
#     to the same cluster (based on ccs_graph_d). Runs on GPU if available.
#     """

#     # --- Device selection ---
#     print("Building graph ",build_graph)
#     device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
#     print("Running on %s" %("GPU" if torch.cuda.is_available() else "CPU"))
#     radius = radius_clustering / 2
#     communities_g = nx.Graph()

#     # Normalize and ensure Torch tensor on chosen device
#     print("normalizing embeddings")
#     t_norm_0 = time.time()
#     embeddings_bincontigs_nz = vamb.cluster._normalize(embeddings_bincontigs)
#     if not torch.is_tensor(embeddings_bincontigs_nz):
#         embeddings_bincontigs_nz = torch.as_tensor(embeddings_bincontigs_nz)
#     embeddings_bincontigs_nz = embeddings_bincontigs_nz.to(device=device, dtype=torch.float32)

#     print("Embeddings normalized in %.2f seconds"%(time.time()-t_norm_0))


#     neighs = [[] for _ in range(len(contignames))]
#     contig_index_lookup = {name: idx for idx, name in enumerate(contignames)}

#     total_clusters = len(ccs_graph_d)
#     fraction_clusters = max(1, round(total_clusters * 0.1))

#     total_contigs = len(contignames)
#     contigs_counter = 0
#     fraction_contigs = max(1, round(len(contignames) * 0.1))
#     t_proces_0 = time.time()
#     with torch.inference_mode():
#         for cluster_n, contigs_in_cluster in enumerate(ccs_graph_d.values()):

#             # Map cluster contig names to global indices (CPU)
#             sorted_contig_idxs = sorted(
#                 contig_index_lookup[c] for c in contigs_in_cluster if c in contig_index_lookup
#             )
#             if len(sorted_contig_idxs) < 2:
#                 contigs_counter += len(sorted_contig_idxs)
#                 if (contigs_counter + 1) > fraction_contigs:
#                     print(f"\t{contigs_counter+1}/{total_contigs} contigs processed {time.time() - t_proces_0}")
#                     fraction_contigs += max(1, round(len(contignames) * 0.1))
#                 if (cluster_n + 1) % fraction_clusters == 0:
#                     print(f"{cluster_n+1}/{total_clusters} clusters processed {time.time() - t_proces_0}")
#                 continue

#             # --- (Minimal change) Build mask on device ---
#             mask_cluster = torch.zeros(len(embeddings_bincontigs_nz), dtype=torch.bool, device=device)
#             idxs_tensor = torch.tensor(sorted_contig_idxs, dtype=torch.long, device=device)
#             mask_cluster[idxs_tensor] = True

#             # Local name list (CPU)
#             sorted_contigs_in_cluster = [contignames[idx] for idx in sorted_contig_idxs]
#             sorted_contigs_in_cluster_to_id = {contig: idx for idx, contig in enumerate(sorted_contigs_in_cluster)}

#             # Iterate contigs in this cluster
#             for i, contig in enumerate(sorted_contigs_in_cluster):
#                 # Compute distances on cluster-filtered embedding array (still your call)
#                 # NOTE: This slices on the device, so calc happens on GPU if available
#                 cluster_distances = vamb.cluster._calc_distances(
#                     embeddings_bincontigs_nz[mask_cluster],
#                     sorted_contigs_in_cluster_to_id[contig]
#                 )

#                 # within-radius mask on device
#                 within_radius_mask = cluster_distances <= radius

#                 # self-distance removal
#                 # i is a Python int for local index; mask is device tensor
#                 if within_radius_mask.numel() > i:
#                     within_radius_mask[i] = False

#                 if not torch.any(within_radius_mask):
#                     continue

#                 # Torch where (no NumPy), then move small indices to CPU once
#                 lat_neigh_local = torch.where(within_radius_mask)[0]
#                 lat_neigh_local_cpu = lat_neigh_local.to("cpu").numpy()  # small transfer

#                 for j, lat_neigh_idx in enumerate(lat_neigh_local_cpu):
#                     # Map local cluster index -> global index (CPU list)
#                     global_neigh_idx = sorted_contig_idxs[lat_neigh_idx]
#                     neighs[contig_index_lookup[contig]].append(global_neigh_idx)

#                     if build_graph:
#                         # Extract distance for this neighbor; small scalar transfer
#                         d = cluster_distances[within_radius_mask][j].item()
#                         communities_g.add_edge(
#                             contig,
#                             sorted_contigs_in_cluster[lat_neigh_idx],
#                             distance=d
#                         )

#             contigs_counter += len(sorted_contig_idxs)
#             if (contigs_counter + 1) > fraction_contigs:
#                 print(f"\t{contigs_counter+1}/{total_contigs} contigs processed in {time.time() - t_proces_0}")
#                 fraction_contigs += max(1, round(len(contignames) * 0.1))

#             if (cluster_n + 1) % fraction_clusters == 0:
#                 print(f"{cluster_n+1}/{total_clusters} clusters processed in {time.time() - t_proces_0}")

#     return (np.array(neighs, dtype=object), communities_g)


# def get_neighbourhoods(neighs_obj, contignames, min_neighs=2):
#     neighbourhoods_g = nx.Graph()
#     for i, neigh_idxs in enumerate(neighs_obj):
#         c = contignames[i]
#         for neigh_idx in neigh_idxs:
#             c_neigh = contignames[neigh_idx]
#             if c_neigh == c:
#                 continue
#             neighbourhoods_g.add_edge(c_neigh, c)

#     hood_cs_d = {
#         i: cc
#         for i, cc in enumerate(nx.connected_components(neighbourhoods_g))
#         if len(cc) >= min_neighs
#     }
#     del neighbourhoods_g
    
#     return hood_cs_d



import time
from collections import defaultdict

import numpy as np
import torch
import networkx as nx

# -----------------------
# Union-Find Structure
# -----------------------
class UnionFind:
    def __init__(self):
        # parent maps integer contig index -> root index
        self.parent = {}

    def find(self, x):
        # lazy init
        if x not in self.parent:
            self.parent[x] = x
        # path compression
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def find_neighbours_optimized(
    embeddings_bincontigs,
    contignames,
    ccs_graph_d,
    radius_clustering=0.2,
    build_graph=False,
    binsplit_separator="C",
    min_neighs=2,
    #deduplicate_neighs=True,
    
):
    """
    Optimized neighbor finder that:
      • computes embedding-radius neighbors per cluster (GPU-accelerated)
      • builds neighs_obj (list of neighbor indices per contig)
      • optionally computes connected components via Union-Find
      • optionally builds a NetworkX graph (if build_graph=True)
      • optionally removes components whose members all have unique sample IDs
        (and cleans neigh lists consistently in both directions)

    Returns (order as you specified):
        neighs_obj, (communities_g if build_graph else None), hood_cs_d
    """

    print("Building graph:", build_graph)
    device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
    print("Running on", "GPU" if torch.cuda.is_available() else "CPU")

    radius = radius_clustering / 2
    communities_g = nx.Graph() if build_graph else None

    # --- Normalize embeddings on chosen device ---
    print("normalizing embeddings")
    t_norm_0 = time.time()
    embeddings_bincontigs_nz = vamb.cluster._normalize(embeddings_bincontigs)
    if not torch.is_tensor(embeddings_bincontigs_nz):
        embeddings_bincontigs_nz = torch.as_tensor(embeddings_bincontigs_nz)
    embeddings_bincontigs_nz = embeddings_bincontigs_nz.to(device=device, dtype=torch.float32)
    print("Embeddings normalized in %.2f seconds" % (time.time() - t_norm_0))

    # neighs_obj storage (global indices)
    neighs = [[] for _ in range(len(contignames))]

    # Union-Find to accumulate components while we discover neighbors
    uf = UnionFind() 

    # Global name -> index lookup
    contig_index_lookup = {name: idx for idx, name in enumerate(contignames)}

    # Progress
    total_clusters = len(ccs_graph_d)
    next_cluster_checkpoint = max(1, round(total_clusters * 0.1))
    total_contigs = len(contignames)
    contigs_counter = 0
    next_contig_checkpoint = max(1, round(total_contigs * 0.1))
    t_proces_0 = time.time()

    with torch.inference_mode():
        for cluster_n, contigs_in_cluster in enumerate(ccs_graph_d.values(), start=1):

            # Map cluster contig names -> global indices
            sorted_contig_idxs = sorted(
                contig_index_lookup[c] for c in contigs_in_cluster if c in contig_index_lookup
            )

            # Small clusters: skip heavy work
            if len(sorted_contig_idxs) < 2:
                contigs_counter += len(sorted_contig_idxs)
                if contigs_counter >= next_contig_checkpoint:
                    print(
                        f"\t{contigs_counter}/{total_contigs} contigs processed {time.time() - t_proces_0:.2f}s"
                    )
                    next_contig_checkpoint += max(1, round(total_contigs * 0.1))
                if cluster_n % next_cluster_checkpoint == 0:
                    print(f"{cluster_n}/{total_clusters} clusters processed {time.time() - t_proces_0:.2f}s")
                continue

            # *** Faster than full-length mask: directly index by the cluster indices ***
            idxs_tensor = torch.tensor(sorted_contig_idxs, dtype=torch.long, device=device)
            cluster_embed = embeddings_bincontigs_nz.index_select(0, idxs_tensor)
            # Local names for graph edges (CPU)
            sorted_contigs_in_cluster = [contignames[idx] for idx in sorted_contig_idxs]

            # Loop through contigs in this cluster
            for i_local, contig_name in enumerate(sorted_contigs_in_cluster):
                # Compute distances on the GPU-sliced cluster embedding matrix
                # Use i_local directly instead of mapping name -> local index
                cluster_distances = vamb.cluster._calc_distances(
                    cluster_embed,  # shape: (#cluster, embedding_dim)
                    i_local         # pivot (row) index inside cluster
                )

                # Within-radius mask
                within_radius_mask = cluster_distances <= radius
                # Remove self; (numel check is unnecessary)
                within_radius_mask[i_local] = False

                # Fast path: no neighbors
                if not torch.any(within_radius_mask):
                    continue

                # Neighbor local indices (device -> CPU small copy)
                neigh_local = torch.where(within_radius_mask)[0]
                neigh_local_cpu = neigh_local.cpu().numpy()

                # Map to global indices and fill outputs
                contig_global_idx = contig_index_lookup[contig_name]

                # Optional: if building graph, fetch distances once efficiently
                if build_graph:
                    # distances at neighbor positions only (no boolean materialization)
                    neigh_dists = cluster_distances.index_select(0, neigh_local).cpu().tolist()

                for j, neigh_local_idx in enumerate(neigh_local_cpu):
                    global_neigh_idx = sorted_contig_idxs[neigh_local_idx]

                    # Store in neighs_obj
                    neighs[contig_global_idx].append(global_neigh_idx)

                    # DSU union (build components incrementally)
                    uf.union(contig_global_idx, global_neigh_idx)

                    # Optional graph edge with distance attribute
                    if build_graph:
                        communities_g.add_edge(
                            contig_name,
                            sorted_contigs_in_cluster[neigh_local_idx],
                            distance=neigh_dists[j]
                        )

            # Progress logging
            contigs_counter += len(sorted_contig_idxs)
            if contigs_counter >= next_contig_checkpoint:
                print(
                    f"\t{contigs_counter}/{total_contigs} contigs processed in {time.time() - t_proces_0:.2f}s"
                )
                next_contig_checkpoint += max(1, round(total_contigs * 0.1))
            if cluster_n % next_cluster_checkpoint == 0:
                print(f"{cluster_n}/{total_clusters} clusters processed in {time.time() - t_proces_0:.2f}s")

    # Build components from DSU if requested
    
    comps = defaultdict(set)
    for idx, name in enumerate(contignames):
        root = uf.find(idx)
        comps[root].add(name)
    hood_cs_d = {
        i: comp
        for i, comp in enumerate(comps.values())
        if len(comp) >= min_neighs
    }

    # Determine which component IDs to drop
    hoods_to_remove = []
    for hood_id, cs in hood_cs_d.items():
        sample_ids = {c.split(binsplit_separator)[0] for c in cs}
        if len(sample_ids) == len(cs):
            hoods_to_remove.append(hood_id)

    if len(hoods_to_remove) > 0:
        for hood in hoods_to_remove:
            cs_to_clear = hood_cs_d[hood]
            for c in cs_to_clear:
                neighs[contig_index_lookup[c]]=[]
            
        # Finally remove those components from the dict
        for hood_id in hoods_to_remove:
            del hood_cs_d[hood_id]

    neighs_obj = np.array(neighs, dtype=object)
    
    return neighs_obj, hood_cs_d, (communities_g if build_graph else None)

if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="From the embeddings, the list of contigs that have been embedded, and the binning contigs, generate a neighbours file, a contigs with neighbours mask, and a dictionary with neighbours and hoods."
    )
    # Add arguments
    parser.add_argument("--embs", help="Embeddings file.")
    parser.add_argument("--contigs_embs", help="Contigs embedded file.")
    parser.add_argument("--neighs_outdir", help="Neighbours outtdir.")
    parser.add_argument("-g", help="Graph files that was used as input for n2v")
    parser.add_argument(
        "--contignames",
        default=None,
        help="File containing the binning contignames, if not provided it will generate neighs with all contigs.",
    )
    parser.add_argument(
        "-r",
        default=[0.1],
        help="Radius within 2 contigs are considered neighbours.",
        nargs="+",
    )

    ## Taking the n2v embeddings from either the assembly graph or the alignment graph. Extract communities of contigs from them.
    # 1. Sort and mask so neighbours are generated according the order defined by the contignames file content
    # 2. Compute cosine distances all against all, only for considered contigs
    # 3. For each contig, define as neighbours the ones that are within a radius
    # 4. For each community of connected contigs, split community if belong to different graph componentnts

    ## Print git commit so we can debug
    commit_hash = get_git_commit(os.path.abspath(__file__))
    print(f"Git commit hash: {commit_hash}")

    ## Parse arguments
    args = parser.parse_args()
    try:
        os.mkdir(args.neighs_outdir)
    except:
        pass

    with open(os.path.join(args.neighs_outdir, "log.txt"), "w") as f:
        f.write("%s\n" % (args))

    ## Load files
    # Load graph file
    t0=time.time()
    with open(args.g, "rb") as pkl_file:
        graph = pickle.load(pkl_file)
    print("Graph loaded with %i nodes and %i edges in %.2f seconds" % (graph.number_of_nodes(), graph.number_of_edges(),time.time()-t0))
    
    # Define which component each contig belongs to
    ccs_graph_d = {i: cc for i, cc in enumerate(nx.connected_components(graph))}

    # Define radiuses
    radiuses = [float(r) for r in args.r]

    # Load contignames
    print("Loading contig names...")
    if args.contignames != None:
        ## Load the names of the contigs that will be used for binning BE CAREFUL SINCE HE CONTIGNAMES USUALLY USED ONLY APPLUES FOR MIN CONTIG LEN 2000
        contignames = np.loadtxt(args.contignames, dtype=object)
    else:
        contignames = np.loadtxt(args.contigs_embs, dtype=object)

    print("len of contignames %i" % (len(contignames)))

    contignames_set = set(contignames)

    ccs_graph_only_binning_contigs_d = {}
    ccs_graph_only_binning_contigs_d_counts = []
    for i, cc in ccs_graph_d.items():
        inter = cc & contignames_set   
        if len(inter) > 1:
            ccs_graph_only_binning_contigs_d[i] = inter
            ccs_graph_only_binning_contigs_d_counts.append(len(inter))
    del ccs_graph_d
    print("Number of clusters with more than 1 contig in the binning set: %i" % len(ccs_graph_only_binning_contigs_d.keys()))
    print("Average (std) contigs per cluster: %.2f (%.2f)" % (np.mean(ccs_graph_only_binning_contigs_d_counts),np.std(ccs_graph_only_binning_contigs_d_counts)))
    
    c_idx_d = {c: i for i, c in enumerate(contignames)}

    # Load the embeddings, and the contigs that are represented in those embeddings
    contigsembs = np.loadtxt(args.contigs_embs, dtype=object)
    embeddings = np.load(args.embs)["arr_0"]
    contig_emb_d = {c: e for c, e in zip(contigsembs, embeddings)}
    

    ## 1. Mask the embeddings and process them so they match the binning contigs
    embeddings_bincontigs = np.zeros((len(contignames), embeddings.shape[1]))
    del embeddings
    embeddings_mask = np.ones(len(contignames), dtype=bool)
    for i, c in enumerate(contignames):
        if c not in contig_emb_d.keys():
            embeddings_mask[i] = False
            continue
        embeddings_bincontigs[c_idx_d[c], :] = contig_emb_d[c]
    del contig_emb_d 

    fraction_embedded_contigs = np.sum(embeddings_mask) / len(contignames)
    print("Fraction of contigs with embeddings %.3f" % (fraction_embedded_contigs))

    ## 2. Compute cosine distances all against all, only for considered AND contigs and contigs taht are embedded
    t0 = time.time()

    ## 3. For each contig, define as neighbours the ones that are within a radius
    embs_d = dict()
    for radius in radiuses:
        print("Finding neighbours within radius %.3f" % radius)
        embs_d[radius] = dict()
        embs_d[radius]["neighs"] = find_neighbours_optimized(
            embeddings_bincontigs,
            contignames,
            ccs_graph_only_binning_contigs_d,
            radius,
            build_graph=False,
            binsplit_separator="C"
        )        

        print("Optimized version finished in %.2f seconds" % (time.time() - t0))
        print("%i contigs with neighs"%(len([n for n in embs_d[radius]["neighs"][0] if len(n) > 0])))
        if args.contignames != None:
            
            neighs_file = os.path.join(
                args.neighs_outdir,
                "neighs_intraonly_rm_object_r_%s.npz" % (str(radius)),
            )
            np.savez_compressed(
                neighs_file,
                embs_d[radius]["neighs"][0],
            )
            print(
                f"Neighs where only intra edges hoods are removed by sample saved in {neighs_file}"
            )

        ## also save cleared hoods where only intra edges hoods are removed
        hoods_clusters_path = os.path.join(
            args.neighs_outdir, "hoods_intraonly_rm_clusters_r_%s.tsv" % (str(radius))
        )

        with open(hoods_clusters_path, "w") as f:
            f.write("clustername\tcontigname\n")
            for hood_i, cs in embs_d[radius]["neighs"][1].items():
                for c in cs:
                    f.write("%s\t%s\n" % (hood_i, c))

        print(
            f"Hoods where only intra edges hoods are removed by sample clusters saved in {hoods_clusters_path}"
        )

    # Save the dictionary to a Pickle file
    # path_dict = os.path.join(embs_dir, "tmp/embs_d_%s.pkl"%(date))
    path_dict = os.path.join("%s/embs_d.pkl" % args.neighs_outdir)

    with open(path_dict, "wb") as pickle_file:
        pickle.dump(embs_d, pickle_file)

    print(f"Dictionary saved to {path_dict}")

