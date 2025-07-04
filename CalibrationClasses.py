from enum import Enum

import scipy.stats as stats

import CalibrationSettings as Sets
import MultiSurvivalModelClasses as SurvivalCls
from SurvivalModelClasses import Cohort
from deampy.calibration import CalibrationRandomSampling, CalibrationMCMCSampling


class CalibrationColIndex(Enum):
    """ indices of columns in the calibration results cvs file  """
    ID = 0          # cohort ID
    W = 1           # likelihood weight
    MORT_PROB = 2   # mortality probability


def log_likelihood(thetas, seed):
    """Compute the log-likelihood of observed data given theta."""

    cohort = Cohort(id=seed, pop_size=Sets.SIM_POP_SIZE, mortality_prob=thetas[0])
    cohort.simulate(n_time_steps=Sets.TIME_STEPS)  # Simulate the cohort

    # get the average survival time for this cohort
    mean = cohort.cohortOutcomes.meanSurvivalTime

    # construct a normal likelihood
    # with mean calculated from the simulated data and standard deviation from the clinical study.
    # evaluate this pdf (probability density function) at the mean reported in the clinical study.
    log_l = stats.norm.logpdf(
        x=Sets.OBS_MEAN,
        loc=mean,
        scale=Sets.OBS_STDEV)

    return log_l

class Calibration:
    def __init__(self):
        """ initializes the calibration object"""
        self.calib = None  # will be set in the sample_posterior method

    def sample_posterior(self, method, n_samples, std_factor=0.1):
        """ sample the posterior distribution of the mortality probability,
        :param method: the sampling method to use, either 'random sampling' or 'mcmc sampling
        :param n_samples: number of samples from the posterior distribution
        :param std_factor: standard deviation factor for the MCMC sampling method
         """

        if method == 'random':
            self.calib = CalibrationRandomSampling(prior_ranges=Sets.PRIOR_RANGE)
            self.calib.run(log_likelihood_func=log_likelihood, num_samples=n_samples)
        elif method == 'mcmc':
            self.calib= CalibrationMCMCSampling(prior_ranges=Sets.PRIOR_RANGE)
            self.calib.run(log_likelihood_func=log_likelihood, std_factor=std_factor, num_samples=n_samples)
        else:
            raise ValueError("Unknown sampling method: {}".format(method))

        self.calib.save_samples(
            file_name="output/samples_{}.csv".format(method),
            parameter_names=['Mortality Probability'])
        if method == 'random':
            self.calib.save_posterior(
                file_name="output/posterior_{}.csv".format(method),
                n_resample=n_samples, parameter_names=['Mortality Probability'], significant_digits=3)


    def plot_posterior(self, method, n_resamples):
        """ plot the posterior distribution of the mortality probability """

        self.calib.read_samples(file_name="output/samples_{}.csv".format(method))
        if method == 'random':
            self.calib.plot_posterior(
                n_resample=n_resamples, figsize=(5, 5), parameter_names=['Mortality Probability'],
                file_name='figs/posterior_plot_{}.png'.format(method))
        elif method == 'mcmc':
            self.calib.plot_posterior(
                n_warmup = int(0.1*n_resamples), figsize=(5, 5), parameter_names=['Mortality Probability'],
                file_name='figs/posterior_plot_{}.png'.format(method))


    # def get_effective_sample_size(self):
    #     """
    #     :returns: the effective sample size
    #     """
    #     return 1 / np.sum(self.normalizedWeights ** 2)


class CalibratedModel:
    """ to run the calibrated survival model """

    def __init__(self, drug_effectiveness_ratio=1):
        """ extracts seeds, mortality probabilities and the associated likelihood from
        the csv file where the calibration results are stored
        :param drug_effectiveness_ratio: effectiveness of the drug
        """

        self.drugEffRatio = drug_effectiveness_ratio

        self.calib = CalibrationRandomSampling(prior_ranges=Sets.PRIOR_RANGE)
        self.calib.read_samples(file_name="output/samples.csv")

        self.multiCohorts = None  # multi-cohort

    def simulate(self, num_of_simulated_cohorts, cohort_size, time_steps):
        """ simulate the specified number of cohorts based on their associated likelihood weight
        :param num_of_simulated_cohorts: number of cohorts to simulate
        :param cohort_size: the population size of cohorts
        :param time_steps: simulation length
        """

        self.calib.resample(n_resample=num_of_simulated_cohorts)

        # simulate the desired number of cohorts
        self.multiCohorts = SurvivalCls.MultiCohort(
            ids=self.calib.resampledSeeds,
            pop_sizes=[cohort_size] * num_of_simulated_cohorts,
            mortality_probs=self.calib.resamples[0] * self.drugEffRatio)

        # simulate all cohorts
        self.multiCohorts.simulate(time_steps)

    def get_mean_survival_time_proj_interval(self, alpha):
        """
        :param alpha: the significance level
        :returns tuple in the form of (mean, [lower, upper]) of projection interval
        """

        mean = self.multiCohorts.multiCohortOutcomes.statMeanSurvivalTime.get_mean()
        proj_interval = self.multiCohorts.multiCohortOutcomes.statMeanSurvivalTime.get_PI(alpha=alpha)

        return mean, proj_interval

