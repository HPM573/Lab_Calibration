import CalibrationClasses as Cls
import CalibrationSettings as Sets
import deampy.plots.histogram as hist
import deampy.plots.sample_paths as path

# initialize a calibrated model
calibrated_model = Cls.CalibratedModel()
# simulate the calibrated model
calibrated_model.simulate(num_of_simulated_cohorts=Sets.NUM_SIM_COHORTS,
                          cohort_size=Sets.SIM_POP_SIZE,
                          time_steps=Sets.TIME_STEPS)

# plot the sample paths
path.plot_sample_paths(
    sample_paths=calibrated_model.multiCohorts.multiCohortOutcomes.survivalCurves,
    title='Survival Curves',
    x_label='Time-Step (Year)',
    y_label='Number Survived',
    transparency=0.25,
    x_range=(0, 100),
    file_name='figs/calibrated/survival_curves.png')

# plot the histogram of mean survival time
hist.plot_histogram(
    data=calibrated_model.multiCohorts.multiCohortOutcomes.meanSurvivalTimes,
    title='Histogram of Mean Survival Time',
    x_label='Mean Survival Time (Year)',
    y_label='Count',
    x_range=[0, 25],
    bin_width=0.5,
    file_name='figs/calibrated/survival_times.png')


# report mean and projection interval
print('Mean survival time and {:.{prec}%} projection interval:'.format(1 - Sets.ALPHA, prec=0),
      calibrated_model.get_mean_survival_time_proj_interval(alpha=Sets.ALPHA))
