import CalibrationClasses as Cls
import CalibrationSettings as Sets
import SimPy.Plots.Histogram as Hist
import SimPy.Plots.SamplePaths as Path


# initialize a calibrated model
calibrated_model = Cls.CalibratedModel(csv_file_name='CalibrationResults.csv')
# simulate the calibrated model
calibrated_model.simulate(num_of_simulated_cohorts=Sets.NUM_SIM_COHORTS,
                          cohort_size=Sets.SIM_POP_SIZE,
                          time_steps=Sets.TIME_STEPS)

# plot the sample paths
Path.plot_sample_paths(
    sample_paths=calibrated_model.multiCohorts.multiCohortOutcomes.survivalCurves,
    title='Survival Curves',
    x_label='Time-Step (Year)',
    y_label='Number Survived',
    transparency=0.5)

# plot the histogram of mean survival time
Hist.plot_histogram(
    data=calibrated_model.multiCohorts.multiCohortOutcomes.meanSurvivalTimes,
    title='Histogram of Mean Survival Time',
    x_label='Mean Survival Time (Year)',
    y_label='Count',
    bin_width=0.25,
    x_range=[2.5, 21.5])

# create the histogram of the resampled mortality probabilities
Hist.plot_histogram(
    data=calibrated_model.resampledMortalityProb,
    title='Histogram of Resampled Mortality Probabilities',
    x_label='Mortality Probability',
    y_label='Counts',
    x_range=[Sets.PRIOR_L, Sets.PRIOR_U])

# Estimate of mortality probability and the posterior interval
print('Estimate of mortality probability ({:.{prec}%} credible interval):'.format(1 - Sets.ALPHA, prec=0),
      calibrated_model.get_mortality_estimate_credible_interval(alpha=Sets.ALPHA))

# report mean and projection interval
print('Mean survival time and {:.{prec}%} projection interval:'.format(1 - Sets.ALPHA, prec=0),
      calibrated_model.get_mean_survival_time_proj_interval(alpha=Sets.ALPHA))
