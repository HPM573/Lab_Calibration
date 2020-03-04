import CalibrationClasses as Cls
import CalibrationSettings as Sets

# create a calibration object
calibration = Cls.Calibration()

# sample the posterior of the mortality probability
calibration.sample_posterior(n_samples=Sets.POST_N)

# effective sample size
txtEff = 'Effective sample size: {:.1f}'.format(calibration.get_effective_sample_size())
print(txtEff)
