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


# # effective sample size
# txtEff = 'Effective sample size: {:.1f}'.format(calibration.get_effective_sample_size())
# print(txtEff)
