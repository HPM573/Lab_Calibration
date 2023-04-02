from enum import Enum

import numpy as np


class CalibrationColIndex(Enum):
    """ indices of columns in the calibration results cvs file  """
    ID = 0          # cohort ID
    W = 1           # likelihood weight
    MORT_PROB = 2   # mortality probability


class Calibration:
    def __init__(self):
        """ initializes the calibration object"""

        self.cohortIDs = []             # IDs of cohorts to simulate
        self.mortalitySamples = []      # values of mortality probability at which the posterior should be sampled
        self.normalizedWeights = []     # normalized likelihood weights (sums to 1)
        self.mortalityResamples = []  # resampled values for constructing posterior estimate and interval

    def sample_posterior(self, n_samples):
        """ sample the posterior distribution of the mortality probability,
         :param n_samples: number of samples from the posterior distribution
         """

        # random number generator

        # cohort ids

        # find values of mortality probability at which the posterior should be evaluated

        # create a multi cohort

        # simulate the multi cohort

        # calculate the likelihood of each simulated cohort

            # get the average survival time for this cohort

            # construct a normal likelihood
            # with mean calculated from the simulated data and standard deviation from the clinical study.
            # evaluate this pdf (probability density function) at the mean reported in the clinical study.

            # store the weight

        # normalize the likelihood weights

        # produce the list to report the results

        # write the calibration result into a csv file

    def get_effective_sample_size(self):
        """
        :returns: the effective sample size
        """
        return 1 / np.sum(self.normalizedWeights ** 2)
