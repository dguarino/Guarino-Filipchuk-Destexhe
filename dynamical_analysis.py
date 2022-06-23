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

# defaults
plt.rcParams['image.cmap'] = 'viridis' # works in old and newer versions of matplotlib

print("... firing statistics")
fr = firinghist(start_time, stop_time, spiketrains, bin_size=frame_duration) #
print("    population firing: {:1.2f}±{:1.2f} sp/frame".format(np.mean(fr),np.std(fr)) )
print("    smoothing")
# Smoothed firing rate: Savitzky-Golay 1D filter
smoothed_fr = savgol_filter(fr, window_length=7, polyorder=3) # as in Filipchuk
smoothed_fr[smoothed_fr<0.] = 0. # filtering can bring the firing below zero
# Baseline: Eilers-Boelens 1D filter
baseline_fr = baseline(smoothed_fr, l=10**8, p=0.01) # as in Filipchuk

fig = plt.figure()
plt.plot(fr, linewidth=0.5, zorder=1)
plt.plot(smoothed_fr, linewidth=0.5, zorder=1)
plt.plot(baseline_fr, linewidth=0.5, zorder=2)
plt.title("mean rate: %.2f (%.2f) sp/frame" % (fr.mean(), fr.std()) )
plt.ylabel('Firing (sp/frame)')
plt.xlabel('Time (frames)')
fig.savefig(exp_path+'/results/firing.svg', transparent=True)
plt.close()
fig.clear()
fig.clf()

print("... generating surrogates to establish population event threshold")
spiketrainsISI = []
cells_cv = []
cells_firing_rate = []
for st in spiketrains:
    cells_firing_rate.append( firinghist_one(st, np.arange(start_time, stop_time, frame_duration) ) ) # cell firing rate
    spiketrainsISI.append( np.diff( st ) )
print("    cells firing rate: {:1.2f}±{:1.2f} sp/s".format(np.mean(cells_firing_rate),np.std(cells_firing_rate)) )

# reshuffle ISIs (100) times
surrogate_fr = []
for isur in range(100):
    # build surrogate rasterplot
    surrogate_spiketrains = []
    for isi in spiketrainsISI:
        random.shuffle(isi) # in-place function
        if len(isi): # init of timing of first spike
            isi[0] += start_time
        surrogate_spiketrains.append( np.cumsum(isi) )
    # compute the population firing histogram for each surrogate timebinned rasterplot
    surrogate_fr.append( firinghist(start_time, stop_time, surrogate_spiketrains, bin_size=frame_duration) )

# instantaneous threshold is the 99% of the surrogate population instantaneous firing rate
event_threshold = np.percentile(np.array(surrogate_fr), 95) + baseline_fr
print("    event size threshold (mean):",np.mean(event_threshold))

print("... find population events in the trial")
# the maximal extrema (beyond threshold) are peaks of population events
peaks = []
for peak in signal.argrelextrema(smoothed_fr, np.greater)[0]:
    if smoothed_fr[peak] < event_threshold[peak]:
        continue # ignore peaks below threshold
    peaks.append(peak)

# the minimal extrema are potential limits of population events
minima = signal.argrelextrema(smoothed_fr, np.less, order=2)[0]
# there can be periods when the firing rate is 0, to be considered minima as well
zerominima = np.where(smoothed_fr == 0)[0]
minima = np.concatenate((minima,zerominima))
minima.sort(kind='mergesort')

# the minimum before and after a peak are taken as start and end times of the population event
events = []
low_toBsaved = []
# each peak (beyond threshold) is part of an event
for peak in peaks:
    event = {'start':0, 'end':0} # init
    # the minima (either below or above threshold) before and after the peak are the real limits of the event
    # start
    for minimum in minima:
        if minimum < peak: # continue until peak
            event['start'] = minimum
            continue
        break # last assigned before peak is taken
    low_toBsaved.append(event['start'])
    # end
    for minimum_idx,minimum in enumerate(minima):
        if minimum <= peak: # continue to peak
            continue
        event['end'] = minima[minimum_idx] # take the one beyond
        break
    low_toBsaved.append(event['end'])
    if (event['start']>=0 and event['end']>0) and (event['end']-event['start']>1):
        events.append(event)
# remove duplicates due to minima above threshold
evti = 0
while evti < len(events):
    evtj = evti + 1
    while evtj < len(events):
        if events[evti] == events[evtj]:
            del events[evtj]
        else:
            evtj += 1
    evti += 1
# removal of unnecessary minima
minima = sorted(low_toBsaved)

# plot everything so far, surrogates, original, threshold, maxima and minima, event
fig, ax = plt.subplots(figsize=(20,5))
x = np.arange(0,len(smoothed_fr))
for event in events:
    ax.axvspan(event['start'], event['end'], alpha=0.1, color='green', zorder=1)
ax.plot(smoothed_fr,linewidth=0.5, color='k', zorder=2)
ax.plot(event_threshold, linewidth=0.5, color='magenta', zorder=3)
plt.plot(baseline_fr, linewidth=0.5, color='orange', zorder=3)
ax.plot(x[minima], smoothed_fr[minima], 'wo', markersize=1, markeredgewidth=0., zorder=4)
ax.plot(x[peaks], smoothed_fr[peaks], 'rs', markersize=.5, zorder=4)
fig.savefig(exp_path+'/results/Events_population_firing.svg', transparent=True)
plt.close()
fig.clear()
fig.clf()

# with not enough events, no point in going further
if len(events)<4:
    print("... not enough events to continue analysis (<4)")
else:
    print("... signatures of population events")
    # produce a cell id signature of each population event event: string of cell ids firing during the event event
    # get the cell ids for each event
    events_signatures = [] # N x M, N events and all M idx for each event (0s and 1s)
    events_spectrums = [] # N x Z, N events and only Z cids for each event
    toBremoved = []
    for event in events:
        # signature
        signature = [0 for i in range(len(spiketrains))] # init
        spectrum = []
        # start end of event index are converted to ms in order to search cell ids being active in that window
        tstart = event['start'] * frame_duration
        tend = event['end'] * frame_duration
        for idx,spiketrain in enumerate(spiketrains):
            # take the idx if there are spikes within event start and end
            spiketrain = np.array(spiketrain)
            # choose idx based on event limits
            if np.any(np.logical_and(spiketrain>=tstart, spiketrain<tend)):
                spectrum.append( ophys_cell_ids[idx] ) # to store the actual cid and not just the index
                signature[idx] = 1 #
        # check that the event signature has more than 1 cell
        if np.count_nonzero(spectrum)>1:
            events_spectrums.append( spectrum )
        if np.count_nonzero(signature)>1:
            events_signatures.append( signature )
        else:
            toBremoved.append(events.index(event))

    # removing events with just one cell
    for index in sorted(toBremoved, reverse=True):
        del events[index]

    events_signatures = np.array(events_signatures)
    events_spectrums = np.array(events_spectrums, dtype=object)

    # Population events number of events per sec
    print("    number of events:",len(events))
    events_sec = len(events_signatures)/stop_time
    print("    number of events per sec:",events_sec)

    # events durations and intervals
    events_durations = []
    events_intervals = []
    last_event = None
    for event in events:
        events_durations.append(event['end']-event['start'])
        if last_event: # only from the second event on
            events_intervals.append(event['start']-last_event['end'])
        last_event = event

    # Population events mean+std Duration
    events_durations_f = np.array(events_durations, dtype=float)
    events_durations_f = events_durations_f*frame_duration
    np.save(exp_path+'/results/events_durations.npy', events_durations_f)
    print("    events duration: %.3f±%.3f" % (np.median(events_durations_f), np.std(events_durations_f)))
    fig = plt.figure()
    plt.yscale('log')
    plt.hist(events_durations_f, bins='auto')
    lims = plt.ylim()
    plt.vlines([np.median(events_durations_f)], ymin=lims[0], ymax=lims[1], linestyles='dashed', colors='k')
    plt.title("median events duration: %.3f±%.3f" % (np.median(events_durations_f), np.std(events_durations_f)) )
    plt.ylabel('event occurrences')
    plt.xlabel('event duration (sec)')
    plt.yscale('log')
    plt.xscale('log')
    fig.savefig(exp_path+'/results/Events_durations.png', transparent=True)
    fig.savefig(exp_path+'/results/Events_durations.svg', transparent=True)
    plt.close()
    fig.clear()
    fig.clf()

    # Population events size (number of cells per event)
    events_size = []
    for esignature in events_signatures:
        events_size.append(len(np.nonzero(esignature)[0]))
    np.save(exp_path+'/results/events_size.npy', events_size)
    print("    events size: %.3f±%.3f" % (np.median(events_size), np.std(events_size)))
    fig = plt.figure()
    plt.hist(events_size, bins='auto')
    lims = plt.ylim()
    plt.vlines([np.median(events_size)], ymin=lims[0], ymax=lims[1], linestyles='dashed', colors='k')
    plt.title("median size: %.3f (%.3f)" % (np.median(events_size), np.std(events_size)) )
    plt.ylabel('occurrences')
    plt.xlabel('number of cells per event')
    plt.yscale('log')
    plt.xscale('log')
    fig.savefig(exp_path+'/results/Events_size.png', transparent=True)
    fig.savefig(exp_path+'/results/Events_size.svg', transparent=True)
    plt.close()
    fig.clear()
    fig.clf()

    # --------------------------------------------------------------------------

    print("... Similarity of events matrix")
    # Raw Pearson's correlation matrix over events signatures (their patterns of cells)
    SimilarityMap = np.corrcoef(events_signatures)
    fig = plt.figure()
    plt.pcolormesh(SimilarityMap)
    cbar = plt.colorbar()
    fig.savefig(exp_path+'/results/Events_CorrMatrix.png', dpi=600, transparent=True)
    plt.close()
    fig.clear()
    fig.clf()

    # --------------------------------------------------------------------------

    # perform clustering linkage by complete cross-correlation of event signatures
    print("... clustering")
    print("    linkage")
    Z = linkage(events_signatures, method='complete', metric='correlation') #
    Z[ Z<0 ] = 0 # for very low correlations, negative values can result
    cut_off = 0.8*max(Z[:,2]) # generic cutoff as in matlab, but we also bootstrap below

    print("    surrogate events signatures for clustering threshold")
    # threshold for cluster significance
    # generate 100 surrogates events_signatures
    # correlate and cluster...
    # there will be clusters, happening just by chance due to the finite number of cells
    # but their internal correlation should not be high
    # the 95% correlation of this random cluster will be the threshold
    # for a cluster to be significantly correlated
    surrogate_reproducibility_list = []
    for csur in range(100):
        surrogate_events_signatures = []
        for evsig in events_signatures:
            surrogate_signature = np.array([0 for i in ophys_cell_indexes]) # init
            surrogate_signature[ np.random.choice(ophys_cell_indexes, size=np.count_nonzero(evsig), replace=False) ] = 1
            surrogate_events_signatures.append(surrogate_signature.tolist())
        # similarity
        surrogate_similaritymap = np.corrcoef(surrogate_events_signatures)
        # clustering
        surrogate_Z = linkage(surrogate_events_signatures, method='complete', metric='correlation') #
        surrogate_Z[ surrogate_Z<0 ] = 0 # for very low correlations, negative values can result
        surrogate_events_assignments = fcluster(surrogate_Z, t=cut_off, criterion='distance')
        # sorting by cluster based on the first element of zip
        surrogate_permutation = [x for _, x in sorted(zip(surrogate_events_assignments, range(len(surrogate_events_signatures))))]
        clustered_surrogate_similaritymap = surrogate_similaritymap[surrogate_permutation] # x
        clustered_surrogate_similaritymap = clustered_surrogate_similaritymap[:,surrogate_permutation] # y
        # cluster reproducibility
        surrogate_events_cluster_sequence = sorted(surrogate_events_assignments) # [ 1 1 1 1 2 2 3 3 3 3 3 ...]
        surrogate_cnt = collections.Counter()
        for word in surrogate_events_cluster_sequence:
            surrogate_cnt[word] += 1
        surrogate_events_cluster_chunks = list(surrogate_cnt.values())
        starti = 0
        for iblock,nblock in enumerate(surrogate_events_cluster_chunks):
            endi = starti+nblock
            # get sub-array
            surrogate_cluster_subarray = np.array(clustered_surrogate_similaritymap[ starti:endi-1, starti:endi-1 ])
            if surrogate_cluster_subarray.size:
                np.fill_diagonal(surrogate_cluster_subarray, 0.0)
                # compute the subarray average
                surrogate_reproducibility_list.append( np.nanmean(surrogate_cluster_subarray) )
            starti = endi

    # statistically significant reproducibility
    cluster_reproducibility_threshold = np.percentile(np.array(surrogate_reproducibility_list), 95)
    print("    cluster reproducibility threshold:",cluster_reproducibility_threshold)
    # number of events in a cluster, even small clusters as long as they pass the reproducibility threshold
    cluster_size_threshold = 2 # minimum requirement
    print("    cluster size threshold:",cluster_size_threshold)

    # clusters
    events_assignments = fcluster(Z, t=cut_off, criterion='distance')
    # print(events_assignments)
    # print(len(events_assignments))
    # [4 4 4 4 2 2 34 4 7 7 7 7 5 5 5 5 4 ... ]

    # count clusters
    c_cnt = collections.Counter()
    for word in events_assignments:
        c_cnt[word] += 1
    nevents_clusters = np.array(list(c_cnt.values()))
    # print("    events/clusters:", nevents_clusters)
    # [ 39  19  11  70  15  74  45  13  10  45 ...]
    nclusters = len(nevents_clusters)
    print("    #clusters:",nclusters)

    # color map of the clustered events
    cmap = mpcm.get_cmap('rainbow')
    cluster_color_array = [mpl.colors.rgb2hex(rgba) for rgba in cmap(np.linspace(0.0, 1.0, nclusters))]
    random.shuffle(cluster_color_array) # to have different nearby colors
    cluster_color_array = np.array(cluster_color_array)
    # print("cluster colors:",len(cluster_color_array))

    threshold_map = nevents_clusters < cluster_size_threshold
    print("    below size threshold:", np.count_nonzero(threshold_map))
    cluster_color_array[threshold_map] = 'gray' # or 'none'

    color_array = []
    for cluidx in events_assignments:
        color_array.append(cluster_color_array[cluidx-1]) # fcluster returns numeric labels from 0 to nclusters-1
    color_array = np.array(color_array)
    # color_array = ['#304030', '#008c85', '#008c85', '#005955', 'gray', 'gray', '#2db3ac', ...]
    events_color_assignments = np.copy(color_array)

    print("    sorting events signatures by cluster")
    # sorting the index of events_signatures based on events_assignments (the first element of zip)
    permutation = [x for _,x in sorted(zip(events_assignments, range(len(events_signatures))))] 
    # permuting all arrays used after
    clustered_signatures = events_signatures[permutation]
    clustered_spectrums = events_spectrums[permutation]
    clustered_event_colors = color_array[permutation]
    events = np.array(events)
    clustered_events = events[permutation]
    np.save(exp_path+'/results/clustered_signatures.npy', clustered_signatures)

    # reordering SimilarityMap
    clustered_SimilarityMap = SimilarityMap[permutation] # x
    clustered_SimilarityMap = clustered_SimilarityMap[:,permutation] # y

    # Pattern reproducibility - Cluster self-similarity
    reproducibility_list = [ 0. for elc in cluster_color_array ]
    core_reproducibility = {elc:1. for elc in cluster_color_array}
    events_cluster_sequence = sorted(events_assignments) # [ 1 1 1 1 2 2 3 3 3 3 3 ...]
    # subarray iteration
    ec_cnt = collections.Counter()
    for word in events_cluster_sequence:
        ec_cnt[word] += 1
    events_cluster_chunks = list(ec_cnt.values())
    starti = 0
    for iblock,nblock in enumerate(events_cluster_chunks):
        endi = starti+nblock
        cluster_subarray = np.array(clustered_SimilarityMap[ starti:endi-1, starti:endi-1 ])
        if cluster_subarray.size:
            np.fill_diagonal(cluster_subarray, 0.0)
            # compute the subarray average
            reproducibility_list[iblock] = np.nanmean(cluster_subarray)
            # overwrite color
            if np.nanmean(cluster_subarray) <= cluster_reproducibility_threshold:
                cluster_color_array[iblock] = "gray"
            else:
                # Stimulus-free method to detect core neurons:
                # within each cluster of events,
                # cores are those participating to more than 95% of cluster events
                core_reproducibility[cluster_color_array[iblock]] = np.percentile(cluster_subarray, 95)
                # core_reproducibility[cluster_color_array[iblock]] = np.percentile(cluster_subarray, 85)
                # core_reproducibility[cluster_color_array[iblock]] = np.percentile(cluster_subarray, 75)
                # core_reproducibility[cluster_color_array[iblock]] = np.percentile(cluster_subarray, 65)
                # core_reproducibility[cluster_color_array[iblock]] = np.percentile(cluster_subarray, 55)
        starti = endi

    # plot all
    fig, ax = plt.subplots()
    plt.pcolormesh(clustered_SimilarityMap)
    # loop over zip(clustered_events,color_array) or over events_assignments
    clcoord = 0
    # for csize,ccolor,reproval in zip(nevents_clusters,clustered_event_colors,reproducibility_list):
    for csize,ccolor,reproval in zip(nevents_clusters,cluster_color_array,reproducibility_list):
        rect = patches.Rectangle((clcoord, clcoord), csize, csize, linewidth=0.5, edgecolor=ccolor, facecolor='none')
        ax.add_patch(rect)
        # ax.text(clcoord+1, clcoord+1, "{:1.2f}".format(reproval), color=ccolor, fontsize=2)
        clcoord += csize
    cbar = plt.colorbar()
    cbar.outline.set_visible(False)
    cbar.set_ticks([]) # remove all ticks
    for spine in plt.gca().spines.values(): # get rid of the frame
        spine.set_visible(False)
    plt.xticks([]) # remove all ticks
    plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
    fig.savefig(exp_path+'/results/Events_CorrClustered.png', transparent=True, dpi=600)
    plt.close()
    fig.clf()

    # Pattern reproducibility by Cluster
    fig = plt.figure()
    for repri, (reprv, reprc) in enumerate(zip(reproducibility_list,cluster_color_array)):
        plt.bar(repri, reprv, 1., facecolor=reprc)
    plt.title('Pattern reproducibility')
    plt.ylabel('Auto-correlation')
    plt.xlabel('Clusters')
    fig.savefig(exp_path+'/results/Pattern_reproducibility.png', transparent=True)
    fig.savefig(exp_path+'/results/Pattern_reproducibility.svg', transparent=True)
    plt.close()
    fig.clear()
    fig.clf()

    print("... finding cluster cores")
    clusters_cores = []
    clusters_cores_by_color = {ecolor:list() for ecolor in clustered_event_colors}
    cluster_color_cores = [[] for clsp in clustered_spectrums]
    currentcl = clustered_event_colors[0]
    cluster_events_list = []
    for cl_idlist, cl_color in zip(clustered_spectrums, clustered_event_colors):
        if cl_color=='gray':
            continue
        # when the color changes, plot the map and reset containers
        if currentcl != cl_color:
            # find common subset of cells in a clusters
            cid_counter = {}
            for event_cids in cluster_events_list:
                for cidc in event_cids:
                    if cidc in cid_counter:
                        cid_counter[cidc] += 1/len(cluster_events_list)
                    else: # create
                        cid_counter[cidc] = 1/len(cluster_events_list)
            cluster_core = []
            for cidkey,cidprob in cid_counter.items():
                # Cores identification, independent of stimuli
                # A cell is considered 'core' of multiple events if it has
                # a frequence of occurrence > core_reproducibility (%) 
                if cidprob >= core_reproducibility[currentcl]:
                    cluster_core.append(cidkey)
                    clusters_cores_by_color[currentcl].append(cidkey)
            clusters_cores.append(cluster_core)
            # reset containers
            cluster_events_list = []
            currentcl = cl_color
        # while the color is the same, append the idx to the current ensemble
        cluster_events_list.append( cl_idlist )

    print("    refining cluster cores")
    # We want to check whether a core is firing unspecifically inside and outside of events
    # For each cluster color, take all its events and group or mask the intervals from the rest of the times.
    # Then for each core of this cluster, avg rate inside events - avg outside events
    np_cells_firing_rate = np.array(cells_firing_rate)
    for caidx,ecolor in enumerate(cluster_color_array):
        if ecolor=='gray':
            continue
        # prepare a row_mask to contain all False for events intervals.
        # It will be used to compute the inside-cluster-events firing rate
        # Its inverse will be used to compute the outsite-cluster-events firing rate
        row_mask = np.array([False for el in range(np_cells_firing_rate.shape[1])])
        # select all events of this cluster
        event_idxs = np.argwhere(events_color_assignments==ecolor).flatten()
        for event_idx in event_idxs:
            # take start and end of the event
            estart = events[event_idx]['start']
            eend = events[event_idx]['end']
            row_mask[ estart:eend ] = True
        # remove cores whose firing rate outside its cluster events is higher than inside
        # because it means they are firing unspecifically
        for coids in clusters_cores_by_color[ecolor]:
            if type(coids)!=type([]): # for single element list retrieved as element
                coids = [coids]
            for coid in coids:
                coidx = ophys_cell_ids.index(coid)
                event_meanrate = np.mean(np_cells_firing_rate[coidx][row_mask])
                outside_meanrate = np.mean(np_cells_firing_rate[coidx])
                meanratediff = event_meanrate-outside_meanrate
                if meanratediff < 0.:
                    if caidx in clusters_cores:
                        if coid in clusters_cores[caidx]:
                            clusters_cores[caidx].remove(coid)

    print("    gathering cores from all clusters ...")
    core_indexes = []
    other_indexes = []
    for dyn_core in clusters_cores:
        core_indexes.extend( [ophys_cell_ids.index(strid) for strid in dyn_core] )
    core_indexes = np.unique(core_indexes)
    print("    # cores:",len(core_indexes))
    # print(core_indexes)
    other_indexes = [i for i in range(len(ophys_cell_ids)) if i not in core_indexes]
    print("    # non-cores:",len(other_indexes))

    print("    plotting single events rasterplots ...")
    source_target_cidx = []
    source_target_color = []
    cores_counts = []
    others_counts = []
    for clidx, (cluster_cids, ecolor) in enumerate(zip(clustered_spectrums, clustered_event_colors)):
        if ecolor=='gray':
            continue
        event_id = events_color_assignments.tolist().index(ecolor) # take the first event of this cluster
        event = events[event_id]
        # take start and end of the ensemble
        estart = (event['start'] * frame_duration )
        eend = (event['end'] * frame_duration )
        # print(ecolor,':', estart,eend)
        # print(cluster_cids)
        # slice cluster cells spiketrains to the start,end interval
        event_spiketrains = []
        event_cidxs = []
        for cid in cluster_cids:
            cidx = ophys_cell_ids.index(cid)
            train = spiketrains[cidx]
            if np.array(train[(train>=estart)*(train<=eend)]).size>0:
                event_spiketrains.append( train[(train>=estart)*(train<=eend)] )
                event_cidxs.append( cidx )
        # print(event_spiketrains)
        if len(event_spiketrains) < np.mean(event_threshold):
            continue

        # sort them based on the fitst element of each and sasve also the last for flow analysis
        sorted_event_cidx = [cidx for _,cidx in sorted(zip(event_spiketrains, event_cidxs), key=lambda ez: ez[0][0])] 
        source_target_cidx.append([sorted_event_cidx[0], sorted_event_cidx[-1]]) # take beginning and end cidx
        source_target_color.append(ecolor)
        # sort them based on the fitst element of each
        event_spiketrains = sorted(event_spiketrains, key=lambda etrain: etrain[0])
        # print(event_spiketrains)

        # print("    plotting spike rasterplot for event from cluster ",ecolor,':', estart,eend)
        fig = plt.figure()
        for row,train in enumerate(event_spiketrains):
            ccol = 'gray'
             # Cores
            if row in core_indexes:
                ccol = 'g'
            plt.scatter( train, [row]*len(train), marker='|', facecolors=ccol, s=150, linewidth=3 )
        plt.ylabel("cell IDs")
        plt.xlabel("time (s)")
        fig.savefig(exp_path+'/results/rasterplot_'+str(clidx)+'.svg', transparent=False, dpi=300)
        plt.close()
        fig.clear()
        fig.clf()

        # plot cells as circles connections as edges
        try:
            fig, ax = plt.subplots()
            for soma_loc, ocid in zip(pyc_ca_soma_loc, ophys_cell_ids):
                ccol = 'lightgray'
                zor = 1
                if ocid in cluster_cids:
                    ccol = 'dimgray'
                    zor = 2
                if ocid in clusters_cores_by_color[ecolor]:
                    ccol = 'g'
                    zor = 3
                ax.scatter( soma_loc[0], soma_loc[1], marker='o', edgecolors=ccol, facecolors=ccol, s=40, zorder=zor )
                # projections from the current cell
                cid_postsyn_list = pyc_ca_syn_df[pyc_ca_syn_df['pre_root_id'] == ocid]['post_root_id'].tolist()
                for ps_id in cid_postsyn_list:
                    ips = ophys_cell_ids.index(ps_id)
                    ax.annotate(
                        "",
                        xy=(soma_loc[0], soma_loc[1]), xycoords='data',
                        xytext=(pyc_ca_soma_loc[ips][0],pyc_ca_soma_loc[ips][1]), textcoords='data',
                        arrowprops=dict(arrowstyle="<-",connectionstyle="arc3",color=ccol),
                        zorder=zor
                    )
            fig.savefig(exp_path+'/results/max_projection_'+str(clidx)+'.svg', transparent=True, dpi=300)
            plt.close()
            fig.clear()
            fig.clf()
        except NameError:
            print("    max projection of the cores per event will not be plotted.")

