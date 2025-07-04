import CalibrationClasses as Cls
import CalibrationSettings as Sets

# create a calibration object
calibration = Cls.Calibration()

# sample the posterior of the mortality probability
# calibration.sample_posterior(method='random', n_samples=Sets.N_SAMPLES)
# calibration.plot_posterior(method='random', n_resamples=Sets.N_SAMPLES)

calibration.sample_posterior(method='mcmc', n_samples=Sets.N_SAMPLES, std_factor=Sets.STD_Factor)
calibration.plot_posterior(method='mcmc', n_resamples=Sets.N_SAMPLES)
calibration.calib.plot_trace(
    figsize=(5, 5), file_name='figs/mcmc_trace.png',
    parameter_names=['Mortality Probability'], moving_ave_window=int(0.1 * Sets.N_SAMPLES))

# # effective sample size
# txtEff = 'Effective sample size: {:.1f}'.format(calibration.get_effective_sample_size())
# print(txtEff)
