# general use imports
import random
random.seed(123456)

import numpy as np
np.random.seed(123456)

import os
import io
import os.path
import shutil
import copy
import pprint
import gc
import json
import glob
import gzip
import urllib
import math
import ast
import sys
import wget
import resource
import collections
import itertools
import more_itertools

from functools import cmp_to_key
from functools import reduce
from itertools import zip_longest # analysis
import pickle
# import pickle5 as pickle


import scipy.io as sio
import scipy.signal as signal
# from scipy.signal import find_peaks
from scipy.signal import savgol_filter
import scipy.cluster.hierarchy as hierarchy
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import to_tree
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist, jaccard, squareform
import scipy.stats as stats
from scipy.stats import ttest_ind
from scipy.stats import shapiro
from scipy.stats import hypergeom
from scipy.stats import entropy
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use('Agg') # to be used when DISPLAY is undefined
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as mpcm
import matplotlib.colors as mpcolors
from matplotlib.colors import rgb2hex
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FormatStrFormatter
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.patches as patches

import pandas as pd
pd.options.display.max_colwidth = 100

import igraph as ig
from igraph import *

# ------------------------------------------------
# classes and functions

# Cross-correlation with maxlag
# from: https://stackoverflow.com/questions/30677241/how-to-limit-cross-correlation-window-width-in-numpy
from numpy.lib.stride_tricks import as_strided
def crosscorrelation(x, y, maxlag, mode='corr'):
    """
    Cross correlation with a maximum number of lags.

    `x` and `y` must be one-dimensional numpy arrays with the same length.

    This computes the same result as
        numpy.correlate(x, y, mode='full')[len(a)-maxlag-1:len(a)+maxlag]

    The return vaue has length 2*maxlag + 1.
    """
    py = np.pad(y.conj(), 2*maxlag, mode='constant')
    T = as_strided(py[2*maxlag:], shape=(2*maxlag+1, len(y) + 2*maxlag), strides=(-py.strides[0], py.strides[0]))
    px = np.pad(x, maxlag, mode='constant')
    if mode == 'dot':       # get lagged dot product
        return T.dot(px)
    elif mode == 'corr':    # gets Pearson correlation
        return (T.dot(px)/px.size - (T.mean(axis=1)*px.mean())) / (np.std(T, axis=1) * np.std(px))

# Finds baseline in the firing rate
# Eilers and H. Boelens 2005
# Memory optimised version: https://stackoverflow.com/questions/29156532/python-baseline-correction-library
# Params:
#   l for smoothness (λ)
#   p for asymmetry
# Both have to be tuned to the data at hand.
# We found that generally is a good choice (for a signal with positive peaks):
#   10^2 ≤ l ≤ 10^9
#   0.001 ≤ p ≤ 0.1
# but exceptions may occur.
# In any case one should vary l on a grid that is approximately linear for log l
from scipy import sparse
from scipy.sparse.linalg import spsolve
def baseline(y, l, p, niter=10):
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    D = l * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    for i in range(niter):
        W.setdiag(w) # Do not create a new matrix, just update diagonal values
        Z = W + D
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z


def firinghist( start, end, spiketrains, bin_size=10 ):
    if len(spiketrains)==0:
        return NaN
    # create bin edges based on start and end of slices and bin size
    bin_edges = np.arange( start, end, bin_size )
    # binning total time, and counting the number of spike times in each bin
    hist = np.zeros( bin_edges.shape[0]-1 )
    for spike_times in spiketrains:
        hist += np.histogram( spike_times, bin_edges )[0]
    return hist # no average over population and bin

def firinghist_one( spiketrain, bin_edges ):
    return np.zeros( bin_edges.shape[0]-1 ) + np.histogram( spiketrain, bin_edges )[0]

def isi( spiketrains, bins ):
    isih = np.zeros(bins)
    for st in spiketrains:
        isih += np.histogram( np.diff( st ), bins )[0]
    return isih

def cv( spiketrains, bins ):
    ii = isi(spiketrains, bins)
    return np.std(ii) / np.mean(ii)

from collections import Counter
from collections import OrderedDict
class OrderedCounter(Counter, OrderedDict):
    'Counter that remembers the order elements are first encountered'
    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, OrderedDict(self))
    def __reduce__(self):
        return self.__class__, (OrderedDict(self),)

# ------------------------------------------------

# MICrONS specific functions

# Dictionary with segment id as key.
# scan : Scan id where the cell’s activity was recorded (stimulus differs by scan).
# trace_raw : Raw calcium trace (scaled version of dF).
# trace : Denoised calcium trace.
# spike : Spike probability inferred from the raw calcium trace.
# stimulus : Stimulus label corresponding to the scan the cell’s activity was recorded.
# Traces have a length of 27100 (Note that the raw data has 27300 frames). First 200 frames (blank stimulus) in the raw data are not used as it was used for the extraction of the traces.

# ------------------------------------------------------------------------------
# From
# https://github.com/AllenInstitute/MicronsBinder/blob/master/notebooks/vignette_analysis/function/lib/data.py
# Get soma center coordinates ([um])
def get_soma_loc(data_dict, seg_id):
    seg_ids = data_dict["segment_id"]
    soma_locs = data_dict["loc"]
    return soma_locs[seg_ids==seg_id][0]

# Get calcium trace
def get_trace(data_dict, seg_id, scan_id, trace_type="spike"):
    seg_ids = data_dict["segment_id"]
    scan_ids = data_dict["scan_id"]
    traces = data_dict[trace_type]
    valid = (seg_ids==seg_id)*(scan_ids==scan_id)
    return traces[np.where(valid)[0][0]]

# Get scan id
def get_scan(data_dict, seg_id):
    seg_ids = data_dict["segment_id"]
    scan_ids = data_dict["scan_id"]
    return scan_ids[seg_ids==seg_id][0]

# Get synapse_conneciton density
def get_density(data_dict, seg_id, dens_type="inconn_dens"):
    seg_ids = data_dict["segment_id"]
    densities = data_dict[dens_type]
    return densities[seg_ids==seg_id][0]

# ------------------------------------------------------------------------------

class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return np.ma.masked_array(np.interp(value, x, y))