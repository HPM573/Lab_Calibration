import CalibrationSettings as Sets
from CalibrationClasses import calibrate_random_sampling, calibrate_mcmc, plot_posterior, plot_mcmc_trace


if __name__ == "__main__":

    # calibrate the model using random sampling
    calibrate_random_sampling(n_samples=Sets.N_SAMPLES)
    # plot the posterior distributions
    plot_posterior(method='random', n_resamples=Sets.N_SAMPLES)

    # calibrate the model using MCMC sampling
    calibrate_mcmc(n_samples=Sets.N_SAMPLES, std_factor=Sets.STD_Factor)
    plot_posterior(method='mcmc', n_resamples=Sets.N_SAMPLES)
    plot_mcmc_trace()

#
# # create a calibration object
# calibration = Calibration()
#
# # use random sampling to sample the posterior distribution
# calibration.sample_posterior(method='random', n_samples=Sets.N_SAMPLES)
# calibration.plot_posterior(method='random', n_resamples=Sets.N_SAMPLES, weighted=True)
#
# # use MCMC sampling to sample the posterior distribution
# calibration.sample_posterior(method='mcmc', n_samples=Sets.N_SAMPLES, std_factor=Sets.STD_Factor)
# calibration.plot_posterior(method='mcmc', n_resamples=Sets.N_SAMPLES)
# calibration.calib.plot_trace(
#     figsize=(5, 5), file_name='figs/mcmc_trace.png',
#     moving_ave_window=int(0.1 * Sets.N_SAMPLES))


# # effective sample size
# txtEff = 'Effective sample size: {:.1f}'.format(calibration.get_effective_sample_size())
# print(txtEff)
