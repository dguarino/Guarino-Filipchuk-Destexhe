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
# or ...

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
plt.scatter( degrees,local_clustering_coefficients, marker='o', facecolor='#111111', s=50, edgecolors='none', alpha=0.5) 
plt.plot(degrees,pfit,c='k')
plt.yscale('log')
plt.xscale('log')
# plt.xlim([10,200])
# plt.ylim([0.01,1])
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

print('... Betweenness centrality')
betweenness_centrality = np.array(dgraph.betweenness(vertices=None, directed=True, cutoff=None, weights=None))
np.save(exp_path+'/results/betweenness_centrality.npy', betweenness_centrality)

print('... Motifs')
# number of motifs by class
motifs = np.array(dgraph.motifs_randesu(size=3, cut_prob=None, callback=None))
# which vertices belong to which class
motif_vertices = {}
def cb_motif(graph,vs_motif,isoc):
    if not isoc in motif_vertices:
        motif_vertices[isoc] = []
    motif_vertices[isoc].append(vs_motif)
dgraph.motifs_randesu(size=3, cut_prob=None, callback=cb_motif)

# generate 100 random networks with the same
# number of vertices and number of edges as dgraph
surrogate_motifs = []
for surrogateg in range(100):
    erg = ig.Graph.Erdos_Renyi(n=dgraph.vcount(), m=dgraph.ecount(), directed=True, loops=False)
    surrogate_motifs.append( erg.motifs_randesu(size=3, cut_prob=None, callback=None) )
    ig.plot(erg, exp_path+'/results/ring_erg.png', layout=erg.layout("circle"), edge_curved=0.2, edge_color='orange', edge_width=0.5, edge_arrow_size=0.1, vertex_size=5, vertex_frame_color='orange', vertex_color='orange', margin=50)
surrogate_motifs = np.percentile(surrogate_motifs, 95, axis=0)
# surrogate motif vertices
surrogate_motif_vertices = {}
def cb_motif(graph,vs_motif,isoc):
    if not isoc in surrogate_motif_vertices:
        surrogate_motif_vertices[isoc] = []
    surrogate_motif_vertices[isoc].append(vs_motif)
# take the last surrogate network (just one for now)
erg.motifs_randesu(size=3, cut_prob=None, callback=cb_motif)

# plot
fig = plt.figure()
plt.bar(range(len(motifs)), motifs, color='k', label='real')
plt.bar(np.array(range(len(motifs)))+0.15, surrogate_motifs, color='orange', label='avg surrogates')
plt.legend()
plt.ylabel('occurrences')
plt.xlabel('motifs types')
fig.savefig(exp_path+'/results/motifs_occurrences.svg', transparent=True)
plt.close()
fig.clear()
fig.clf()
# plot
motifsratio = motifs/surrogate_motifs
fig = plt.figure()
plt.bar(range(len(motifs)), motifsratio, color='k')
plt.ylabel('count relative to random')
plt.xlabel('motifs types')
fig.savefig(exp_path+'/results/motifs_ratio.png', transparent=True)
fig.savefig(exp_path+'/results/motifs_ratio.svg', transparent=True)
plt.close()
fig.clear()
fig.clf()
