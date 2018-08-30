from builtins import range

import numpy as np
import os

from collections import OrderedDict
from scipy.sparse import csr_matrix, dia_matrix, coo_matrix, triu, tril
from scipy.sparse import vstack as sparse_vstack
from scipy.sparse import hstack as sparse_hstack
#import tables
from intervaltree import IntervalTree, Interval
from hictracks.utilities import toBytes
from hictracks.utilities import toString
from hictracks.utilities import check_chrom_str_bytes

import gzip

import cooler

import logging
log = logging.getLogger(__name__)

import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.filterwarnings(action="ignore", message="numpy.dtype size changed")
warnings.filterwarnings(action="ignore", message="numpy.ndarray size changed")
warnings.simplefilter(action='ignore', category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=ImportWarning)
#warnings.simplefilter(action='ignore', category=tables.exceptions.FlavorWarning)


# try to import pandas if exists
try:
    import pandas as pd
    pandas = True
except ImportError:
    pandas = False


class hiCMatrix:
    """
    Class to handle Hi-C matrices
    contains routines to get intrachromosomal distances
    get sub matrices by chrname.
    """

    def __init__(self, matrixFile):
        self.cooler_file = cooler.Cooler(matrixFile)        

    def getBinSize(self):
        """
        estimates the bin size. In case the bin size
        is not equal for all bins (maybe except for the
        bin at the en of the chromosomes) a warning is issued.
        In case of uneven bins, the median is returned.
        """
        if self.bin_size is None:
            chrom, start, end, extra = zip(*self.cut_intervals)
            median = int(np.median(np.diff(start)))
            diff = np.array(end) - np.array(start)

            # check if the bin size is
            # homogeneous
            if len(np.flatnonzero(diff != median)) > (len(diff) * 0.01):
                self.bin_size_homogeneous = False
                if self.non_homogeneous_warning_already_printed is False:
                    log.warning('Bin size is not homogeneous. \
                                      Median {}\n'.format(median))
                    self.non_homogeneous_warning_already_printed = True
            self.bin_size = median
        return self.bin_size

    def get_sub_matrix(self, chrom_region, start_region, end_region):
        return  self.cooler_file.matrix(balance = True).fetch((chrom_region, start_region, end_region))

    def get_start_pos(self, chrom_region, start_region, end_region):
        return self.cooler_file.bins().fetch((chrom_region, start_region, end_region))['start']
        