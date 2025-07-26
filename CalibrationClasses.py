import scipy.stats as stats

import CalibrationSettings as Sets
import MultiSurvivalModelClasses as SurvivalCls
from SurvivalModelClasses import Cohort
from deampy.calibration import CalibrationRandomSampling, CalibrationMCMCSampling


def log_likelihood_func(thetas, seed):
    """Compute the log-likelihood of observed data given theta."""

    cohort = Cohort(id=seed, pop_size=Sets.SIM_POP_SIZE, mortality_prob=thetas[0])
    cohort.simulate(n_time_steps=Sets.TIME_STEPS, seed=seed)  # Simulate the cohort

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


def calibrate_random_sampling(n_samples):
    """Run random sampling calibration with the specified prior ranges and log-likelihood function."""
    calib = CalibrationRandomSampling(prior_ranges=Sets.PRIOR_RANGE)
    calib.run(log_likelihood_func=log_likelihood_func, num_samples=n_samples)
    calib.save_samples(file_name="output/samples_random.csv")


def calibrate_mcmc(n_samples, std_factor):
    """Run MCMC sampling calibration with the specified prior ranges and log-likelihood function."""
    calib = CalibrationMCMCSampling(prior_ranges=Sets.PRIOR_RANGE)
    calib.run(log_likelihood_func=log_likelihood_func, std_factor=std_factor, num_samples=n_samples)
    calib.save_samples(file_name="output/samples_mcmc.csv")


def plot_posterior(method, n_resamples, weighted=True):
    """ plot the posterior distribution of the mortality probability """

    if method == 'random':
        calib = CalibrationRandomSampling(prior_ranges=Sets.PRIOR_RANGE)
        calib.read_samples(file_name="output/samples_random.csv")
        calib.plot_posterior(
            n_resample=n_resamples, weighted=weighted, figsize=(5, 5),
            file_name='figs/calibration/{}_postr.png'.format(method))

    elif method == 'mcmc':
        calib = CalibrationMCMCSampling(prior_ranges=Sets.PRIOR_RANGE)
        calib.read_samples(file_name="output/samples_mcmc.csv")
        calib.plot_posterior(
            n_warmup = int(0.1*n_resamples), figsize=(5, 5), parameter_names=['Mortality Probability'],
            file_name='figs/calibration/{}_postr.png'.format(method))

    else:
        raise ValueError("Unknown calibration method: {}".format(method))


def plot_mcmc_trace():
    """ plot the trace of the MCMC sampling """
    calib = CalibrationMCMCSampling(prior_ranges=Sets.PRIOR_RANGE)
    calib.read_samples(file_name="output/samples_mcmc.csv")
    calib.plot_trace(
        figsize=(5, 5), file_name='figs/calibration/mcmc_trace.png',
        moving_ave_window=int(0.1 * Sets.N_SAMPLES))


class CalibratedModel:
    """ to run the calibrated survival model """

    def __init__(self, calib_method, drug_effectiveness_ratio=1):
        """ extracts seeds, mortality probabilities and the associated likelihood from
        the csv file where the calibration results are stored
        :param drug_effectiveness_ratio: effectiveness of the drug
        """

        self.drugEffRatio = drug_effectiveness_ratio

        if calib_method == 'random':
            self.calib = CalibrationRandomSampling(prior_ranges=Sets.PRIOR_RANGE)
        elif calib_method == 'mcmc':
            self.calib = CalibrationMCMCSampling(prior_ranges=Sets.PRIOR_RANGE)
        else:
            raise ValueError("Unknown calibration method: {}".format(calib_method))

        self.calib.read_samples(file_name="output/samples_{}.csv".format(calib_method))

        self.multiCohorts = None  # multi-cohort

    def simulate(self, num_of_simulated_cohorts, cohort_size, time_steps):
        """ simulate the specified number of cohorts based on their associated likelihood weight
        :param num_of_simulated_cohorts: number of cohorts to simulate
        :param cohort_size: the population size of cohorts
        :param time_steps: simulation length
        """

        if isinstance(self.calib, CalibrationRandomSampling):
            # resample the seeds and mortality probabilities
            self.calib.resample(n_resample=num_of_simulated_cohorts)

            seeds = self.calib.resampledSeeds

            # simulate the desired number of cohorts
            self.multiCohorts = SurvivalCls.MultiCohort(
                ids=range(num_of_simulated_cohorts),
                pop_sizes=[cohort_size] * num_of_simulated_cohorts,
                mortality_probs=self.calib.resamples['Mortality Probability'] * self.drugEffRatio)

        elif isinstance(self.calib, CalibrationMCMCSampling):

            seeds = self.calib.seeds[-num_of_simulated_cohorts:]

            # simulate the desired number of cohorts
            self.multiCohorts = SurvivalCls.MultiCohort(
                ids=range(num_of_simulated_cohorts),
                pop_sizes=[cohort_size] * num_of_simulated_cohorts,
                mortality_probs=self.calib.samples['Mortality Probability'][-num_of_simulated_cohorts:] * self.drugEffRatio)
        else:
            raise ValueError("Unknown calibration method: {}".format(self.calib))

        # simulate all cohorts
        self.multiCohorts.simulate(time_steps, seeds=seeds)

    def get_mean_survival_time_proj_interval(self, alpha):
        """
        :param alpha: the significance level
        :returns tuple in the form of (mean, [lower, upper]) of projection interval
        """

        mean = self.multiCohorts.multiCohortOutcomes.statMeanSurvivalTime.get_mean()
        proj_interval = self.multiCohorts.multiCohortOutcomes.statMeanSurvivalTime.get_PI(alpha=alpha)

        return mean, proj_interval

