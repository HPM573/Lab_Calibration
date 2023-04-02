from enum import Enum

import deampy.in_out_functions as IO
import numpy as np
import scipy.stats as stats
from numpy.random import RandomState

import CalibrationSettings as Sets
import MultiSurvivalModelClasses as SurvivalCls


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
        rng = RandomState(1)

        # cohort ids
        self.cohortIDs = range(n_samples)

        # find values of mortality probability at which the posterior should be evaluated
        self.mortalitySamples = rng.uniform(
            low=Sets.PRIOR_L,
            high=Sets.PRIOR_U,
            size=Sets.PRIOR_N)

        # create a multi cohort
        multi_cohort = SurvivalCls.MultiCohort(
            ids=self.cohortIDs,
            mortality_probs=self.mortalitySamples,
            pop_sizes=[Sets.SIM_POP_SIZE] * Sets.PRIOR_N
        )

        # simulate the multi cohort
        multi_cohort.simulate(n_time_steps=Sets.TIME_STEPS)

        # calculate the likelihood of each simulated cohort
        weights = []
        for cohort_id in self.cohortIDs:

            # get the average survival time for this cohort
            mean = multi_cohort.multiCohortOutcomes.meanSurvivalTimes[cohort_id]

            # construct a normal likelihood
            # with mean calculated from the simulated data and standard deviation from the clinical study.
            # evaluate this pdf (probability density function) at the mean reported in the clinical study.
            weight = stats.norm.pdf(
                x=Sets.OBS_MEAN,
                loc=mean,
                scale=Sets.OBS_STDEV)

            # store the weight
            weights.append(weight)

        # normalize the likelihood weights
        sum_weights = sum(weights)
        self.normalizedWeights = np.divide(weights, sum_weights)

        # produce the list to report the results
        csv_rows = \
            [['Cohort ID', 'Likelihood Weights', 'Mortality Prob']]  # list containing the calibration results
        for i in range(len(self.mortalitySamples)):
            csv_rows.append(
                [self.cohortIDs[i], self.normalizedWeights[i], self.mortalitySamples[i]])

        # write the calibration result into a csv file
        IO.write_csv(
            file_name='CalibrationResults.csv',
            rows=csv_rows)

    def get_effective_sample_size(self):
        """
        :returns: the effective sample size
        """
        return 1 / np.sum(self.normalizedWeights ** 2)


class CalibratedModel:
    """ to run the calibrated survival model """

    def __init__(self, csv_file_name):
        """ extracts seeds, mortality probabilities and the associated likelihood from
        the csv file where the calibration results are stored
        :param csv_file_name: name of the csv file where the calibrated results are stored
        """

        # read the columns of the csv files containing the calibration results

        # store likelihood weights, cohort IDs and sampled mortality probabilities

        self.multiCohorts = None  # multi-cohort

    def simulate(self, num_of_simulated_cohorts, cohort_size, time_steps):
        """ simulate the specified number of cohorts based on their associated likelihood weight
        :param num_of_simulated_cohorts: number of cohorts to simulate
        :param cohort_size: the population size of cohorts
        :param time_steps: simulation length
        """
        # resample cohort IDs and mortality probabilities based on their likelihood weights
        # sample (with replacement) from indices [0, 1, 2, ..., number of weights] based on the likelihood weights

        # use the sampled indices to populate the list of cohort IDs and mortality probabilities

        # simulate the desired number of cohorts

        # simulate all cohorts

    def get_mean_survival_time_proj_interval(self, alpha):
        """
        :param alpha: the significance level
        :returns tuple in the form of (mean, [lower, upper]) of projection interval
        """

        return mean, proj_interval

    def get_mortality_estimate_credible_interval(self, alpha):
        """
        :param alpha: the significance level
        :returns tuple (mean, [lower, upper]) of the posterior distribution"""

        # calculate the credible interval


        return estimate, credible_interval
