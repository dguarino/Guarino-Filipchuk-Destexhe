#!/usr/bin/python
"""
Copyright (c) 2022, Domenico GUARINO
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Universite Paris Saclay nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GUARINO BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

random.seed(123456)
np.random.seed(123456)

#--------------------------------------------------------------------------

# if an adjacency matrix is already available, e.g. from electron microscopy, it will be loaded
# in case it is not available, e.g. from calcium imaging, it will be estimated
# from spiking 1-lag correlations, as described in Sadovsky and MacLean 2013.
print("... adjacency matrix")
if os.path.exists(exp_path+'/results/adjacency_matrix.npy'):
    adjacency_matrix = np.load(exp_path+'/results/adjacency_matrix.npy', allow_pickle=True)
    print("... loaded")
else:
    # make binary spiketrains
    print("    creating empty binary spike lists")
    binary_spiketrains = np.zeros( (len(spiketrains),len(time)+2) )
    # print(binary_spiketrains.shape)
    print("    filling binary spike lists")
    # iterate over spiketrains assigning 1 to the binary_spiketrains at the corresponding position
    for row,train in enumerate(spiketrains):
        tidxs = np.trunc(np.array(train)/frame_duration).astype(int)
        tidxs[tidxs>len(time)] = len(time) # double check. The frame_duration is periodic but we represent it to the fourth decimal position.
        binary_spiketrains[row][tidxs] = 1
        # binary_spiketrains[row][np.trunc(train/frame_duration).astype(int)] = 1
    # print(binary_spiketrains)
    print("    composing the adjacency matrix")
    adjacency_matrix = []
    for irow,bsti in enumerate(binary_spiketrains):
        row_xcorr = []
        for jrow,bstj in enumerate(binary_spiketrains):
            if irow==jrow:
                row_xcorr.append(0.0) # no self connections
                continue
            row_xcorr.append(crosscorrelation(bsti, bstj, maxlag=1, mode='corr')[2]) # correlation at lag 1
        adjacency_matrix.append(row_xcorr)
    adjacency_matrix = np.array(adjacency_matrix)
    np.save(exp_path+'/results/adjacency_matrix.npy', adjacency_matrix)
    # plot
    fig = plt.figure()
    plt.pcolormesh(adjacency_matrix)
    cbar = plt.colorbar()
    fig.savefig(exp_path+'/results/adjacency_matrix.png', transparent=True)
    plt.close()
    fig.clear()
    fig.clf()

# The following line only works for inferred connections based on 1-lag correlations.
# The adjaceny matrix from electron microscopy comes with either 0 or 1 connections,
# therefore this line has no effect on EM adjacency matrices.
# To ensure sparseness of the adjacency matrix, SadovskyMacLean2013 discard weak correlations (<0.4)
adjacency_matrix[ adjacency_matrix <= adjacency_matrix.max()*0.4 ] = 0.0
# EM-Ca identified cells (112)
dgraph = ig.Graph.Weighted_Adjacency(adjacency_matrix)

#### or ...

# # Full network
# # igraph import from pandas works better with 1M edges using directly IDs
# # get all root_ids
# root_ids = pd.concat([syn_df["post_root_id"], syn_df["pre_root_id"]], ignore_index=True).drop_duplicates()
# root_idxs = list(range( len(root_ids) )) # all cells
# map_root_ids = dict(zip(root_ids,root_idxs))
# source_s = syn_df["pre_root_id"].map(map_root_ids)
# target_s = syn_df["post_root_id"].map(map_root_ids)
# syn_edges_df = pd.DataFrame({'source':source_s.astype(int),'target':target_s.astype(int)})
# vertices_df = pd.DataFrame(root_idxs, columns=['ID'])

# proofread network
# get all root_ids
root_ids = pd.concat([syn_spines_df["post_root_id"], syn_spines_df["pre_root_id"]], ignore_index=True).drop_duplicates()
root_idxs = list(range( len(root_ids) )) # all cells
map_root_ids = dict(zip(root_ids,root_idxs))
source_s = syn_spines_df["pre_root_id"].map(map_root_ids)
target_s = syn_spines_df["post_root_id"].map(map_root_ids)
syn_edges_df = pd.DataFrame({'source':source_s.astype(int),'target':target_s.astype(int)})
vertices_df = pd.DataFrame(root_idxs, columns=['ID'])

print(len(map_root_ids))
# print(source_s)
# print(len(source_s))
# print(len(target_s))
print(syn_edges_df)
print(vertices_df)
dgraph = ig.Graph.DataFrame(edges=syn_edges_df, directed=True, vertices=vertices_df, use_vids=True)

# ig.plot(dgraph, exp_path+'/results/ring.png', layout=dgraph.layout("circle"), edge_curved=0.2, edge_color='#000', edge_width=0.5, edge_arrow_size=0.1, vertex_size=5, vertex_color='#000', margin=50)

print("    number of vertices:", dgraph.vcount())

print('... Network nodes degrees')
degrees = np.array(dgraph.degree())
print("    ", np.count_nonzero(degrees))
np.save(exp_path+'/results/degrees.npy', degrees)

print("... Degree distributions")
# https://igraph.org/python/api/latest/igraph._igraph.GraphBase.html#degree
degdist = dgraph.degree_distribution(bin_width=5)
degree_counts = [bi[2] for bi in degdist.bins()]
np.save(exp_path+'/results/degree_counts.npy', degree_counts)
fig = plt.figure()
plt.plot(range(len(degree_counts)), degree_counts, linewidth=3.0)
plt.ylabel('Number of vertices')
plt.xlabel('Degree')
plt.xscale('log')
plt.yscale('log')
plt.savefig(exp_path+'/results/degree_distribution.png', transparent=True, dpi=300)
plt.close()
fig.clf()

# Clustering Coefficient of only excitatory cells
print('... Local Clustering Coefficient')
# undirected only
local_clustering_coefficients = np.array(dgraph.transitivity_local_undirected(vertices=None, mode="zero"))
print("    min", np.min(local_clustering_coefficients))
print("    mean", np.mean(local_clustering_coefficients))
print("    max", np.max(local_clustering_coefficients))
np.save(exp_path+'/results/local_clustering_coefficients.npy', local_clustering_coefficients)

# covariance matrix and eigenvalues/vectors
Cov = np.cov(degrees,local_clustering_coefficients)
print("    Covariance matrix degrees,local_clustering_coefficients:")
print(Cov)
w, v = np.linalg.eig(Cov)
print("    eigenvalues:")
print(w)
print("    eigenvectors:")
print(v)

# power-law fit
def powerlaw(x, a, b):
    return a * (x**-b)
# fitparams, _ = curve_fit(powerlaw, degrees, local_clustering_coefficients, p0=np.asarray([1000.,2.]))
# print("fit:",fitparams)
# increasing a pushes the curve up, increasing b tilt clockwise
# paramsfit = [2, 1.5] # 1M EM-only
paramsfit = [2, 1.05] # 334 EM-only all proofread
# paramsfit = [1, 0.5] # 112 Ca/EM
pfit = powerlaw(degrees, *paramsfit)

# Coefficient of Determination
# In the best case, the modeled values exactly match the observed values, which results in ss_res=0 and r2=1.
# A baseline model, which always predicts np.mean(LCC), will have r2=0.
# Models that have worse predictions than this baseline will have a negative r2.
# residual sum of squares
ss_res = np.sum((local_clustering_coefficients - pfit) ** 2)
# total sum of squares
ss_tot = np.sum((local_clustering_coefficients - np.mean(local_clustering_coefficients)) ** 2)
# r-squared goodeness-of-fit
r2 = 1 - (ss_res / ss_tot)

# Hierarchical modularity as in SadovskyMacLean2013
# demonstrated by a linear log-log covariance relationship between node degree and node local clustering coefficient.
fig = plt.figure()
summer = mpcm.summer
# for deg,ccoef in zip(degrees,local_clustering_coefficients):
#     plt.scatter( deg, ccoef, marker='o', facecolor='#AAD400', s=10, edgecolors='none', alpha=0.25) # 22um
plt.scatter( degrees,local_clustering_coefficients, marker='o', facecolor='#111111', s=50, edgecolors='none', alpha=0.5) # 22um
plt.plot(degrees,pfit,c='k')
plt.yscale('log')
plt.xscale('log')
# plt.xlim([10,200])
# plt.ylim([0.1,1])
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.ylabel('LCC')
plt.xlabel('degree')
# ax.set_xticklabels([])
# ax.set_yticklabels([])
plt.tick_params(axis='both', bottom='on', top='on', left='off', right='off')
plt.title("eigv: {:.2f}, {:.4f} - fit x^: {:.2f} - GoF R2={:.2f}".format(w[0], w[1], paramsfit[1], r2))
plt.tight_layout()
fig.savefig(exp_path+'/results/hierarchical_modularity.png', transparent=True, dpi=900)
fig.savefig(exp_path+'/results/hierarchical_modularity.svg', transparent=True)
plt.close()
fig.clf()

print("... MICrONS bow-tie analysis")
# Local bow-ties analysis as in FujitaKichikawaFujiwaraSoumaIyetomi2019
# identify communities based on (multiple trials) random walks as flow
# igraph gets stuck, using the same (and more powerful with teleportation) algorithm
from infomap import Infomap # with teleportation to ensure no local solution
im = Infomap(no_self_links=True, flow_model="directed", seed=2**32-1, core_loop_limit=10, prefer_modular_solution=True, inner_parallelization=True, num_trials=10)
im.add_networkx_graph( dgraph.to_networkx() ) # infomap accepts only networkx format
print("    starting infomap analysis")
im.run()
print(f"    found {im.num_top_modules} modules with codelength: {im.codelength:.4f}  entropy: {im.entropy_rate:.4f}")
previous_id = 1
communities_tot = []
communities_lens = []
community = []
for node_id, module_id in sorted(im.modules, key=lambda x: x[1]):
    if module_id>previous_id: # simple module handling
        community_graph = dgraph.subgraph(community) # community contains the indexes in dgraph
        imcommunity = Infomap(no_self_links=True, flow_model="directed", seed=2**32-1, core_loop_limit=10, prefer_modular_solution=True, silent=True, num_trials=10)
        imcommunity.add_networkx_graph( community_graph.to_networkx() )
        imcommunity.run()
        # for communitynode_id, communitymodule_id in sorted(imcommunity.modules, key=lambda x: x[1]):
        #     print(communitynode_id, communitymodule_id)
        if imcommunity.num_non_trivial_top_modules > 2:
            communities_lens.append(len(community))
        communities_tot.append(len(community))
        # simple module handling
        previous_id=module_id
        community = []
    # print(node_id, module_id)
    community.append(node_id)
print("    bow-tie score:", len(communities_lens)/len(communities_tot))
print("    communities lens:",stats.describe(communities_lens))

print("... local bow-tie analysis")
# Local bow-ties analysis as in FujitaKichikawaFujiwaraSoumaIyetomi2019
# identify communities based on (multiple trials) random walks as flow
# there should not be problems for igraph:
# https://stackoverflow.com/questions/70126260/maximum-amount-of-data-that-r-igraph-package-can-handle
communities = dgraph.community_infomap(trials=100)
print("    communities:",len(communities))
structural_cores = []
communities_lens = []
dgraph_btlabels = np.array( ["#999"] * len(dgraph.vs) )
for icomm,community in enumerate(communities):
    print("    ",icomm,len(community))
    community_btlabels = np.array( ["#999"] * len(community) )
    # create a subgraph to find the bow-tie core
    community_graph = dgraph.subgraph(community) # community contains the indexes in dgraph
    # 1. Find the largest component of this subgraph based on flow
    sorted_subgroups = sorted(community_graph.community_infomap(trials=20), key=len, reverse=True)
    if len(sorted_subgroups)<3:
        continue
    largest = sorted_subgroups[0]
    # get nodes as dgraph vertex indexes
    largest_indexes = np.array(community)[largest].tolist()
    # check that is not alone, in which case we do not consider for bow-tie
    if len(largest)==len(community):
        continue
    # otherwise, let's analyze whether there is a bow-tie structure
    # 2. For each node not in the first, check whether it can reach the first.
    otherS_ids = [sid for sid in range(len(community)) if sid not in largest]
    incomponent = []
    for notcore in community_graph.vs.select(otherS_ids):
        inpaths = community_graph.get_all_simple_paths(notcore, to=largest, cutoff=1, mode='out')
        if len(inpaths)>0:
            incomponent.append(notcore.index)
    # 3. For each node not in the largest nor the in-component
    complement = [sid for sid in otherS_ids if sid not in incomponent]
    # for ocid in complement:
    #     community_btlabels[ocid] = "#11F"
    #     dgraph_btlabels[np.array(community)[ocid]] = "#11F"
    outcomponent = []
    for notcore in community_graph.vs.select(complement):
        outpaths = community_graph.get_all_simple_paths(notcore, to=largest, cutoff=1, mode='in')
        if len(outpaths)>0:
            outcomponent.append(notcore.index)
    # if there is no outcomponent is not a bow-tie
    if len(outcomponent)==0:
        continue
    # in case the bow-tie has been identified, color its nodes and append its community to the count
    communities_lens.append(len(community))
    for compidx in largest:
        dgraph_btlabels[np.array(community)[compidx]] = "#F11"
        community_btlabels[compidx] = "#F11"
    for compidx in incomponent:
        dgraph_btlabels[np.array(community)[compidx]] = "#1F1"
        community_btlabels[compidx] = "#1F1"
    for compidx in outcomponent:
        dgraph_btlabels[np.array(community)[compidx]] = "#11F"
        community_btlabels[compidx] = "#11F"
    # get all dgraph community cores to later check their overlap with dynamical cores
    structural_cores.append( largest_indexes )
    # local community layout
    community_graph.vs["color"] = community_btlabels
    community_graph.es["color"] = [community_graph.vs[edge.source]["color"]+"4" for edge in community_graph.es]
    ig.plot(community_graph, exp_path+'/results/fr_community_'+str(icomm)+'.png', layout=community_graph.layout("fr"), vertex_size=15, vertex_frame_width=0, edge_arrow_size=0.1, margin=50)
print("    communities lens:",stats.describe(communities_lens))
# dgraph Layout
matplotlib.rcParams['figure.dpi'] = 900
# cairo dpi
dgraph.vs["color"] = dgraph_btlabels
dgraph.es["color"] = [dgraph.vs[edge.source]["color"]+"4" for edge in dgraph.es]
ig.plot(dgraph, exp_path+'/results/fr.png', layout=dgraph.layout("mds"), vertex_size=5, vertex_frame_width=0, edge_arrow_size=0.1, margin=50)
