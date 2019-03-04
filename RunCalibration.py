import CalibrationClasses as Cls
import CalibrationSettings as CalibSets
import SimPy.FigureSupport as Fig

# create a calibration object
calibration = Cls.Calibration()

# sample the posterior of the mortality probability
calibration.sample_posterior(n_samples=CalibSets.POST_N)

# create the histogram of the resampled mortality probabilities
Fig.graph_histogram(
    data=calibration.mortalityResamples,
    title='Histogram of Resampled Mortality Probabilities',
    x_label='Mortality Probability',
    y_label='Counts',
    x_range=[CalibSets.POST_L, CalibSets.POST_U])

# Estimate of mortality probability and the posterior interval
print('Estimate of mortality probability ({:.{prec}%} credible interval):'.format(1-CalibSets.ALPHA, prec=0),
      calibration.get_mortality_estimate_credible_interval(alpha=CalibSets.ALPHA))

# effective sample size
txtEff = 'Effective sample size: {:.1f}'.format(calibration.get_effective_sample_size())
print(txtEff)
