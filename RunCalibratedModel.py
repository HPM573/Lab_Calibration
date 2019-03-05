import CalibrationClasses as Cls
import CalibrationSettings as CalibSets
import SimPy.FigureSupport as Fig
import SimPy.SamplePathClasses as Path


# initialize a calibrated model
calibrated_model = Cls.CalibratedModel(csv_file_name='CalibrationResults.csv')
# simulate the calibrated model
calibrated_model.simulate(num_of_simulated_cohorts=CalibSets.NUM_SIM_COHORTS,
                          cohort_size=CalibSets.SIM_POP_SIZE,
                          time_steps=CalibSets.TIME_STEPS)

# plot the sample paths
Path.graph_sample_paths(
    sample_paths=calibrated_model.multiCohorts.multiCohortOutcomes.survivalCurves,
    title='Survival Curves',
    x_label='Time-Step (Year)',
    y_label='Number Survived',
    transparency=0.5)

# plot the histogram of mean survival time
Fig.graph_histogram(
    data=calibrated_model.multiCohorts.multiCohortOutcomes.meanSurvivalTimes,
    title='Histogram of Mean Survival Time',
    x_label='Mean Survival Time (Year)',
    y_label='Count',
    bin_width=0.25,
    x_range=[2.5, 21.5])

# report mean and projection interval
print('Mean survival time and {:.{prec}%} projection interval:'.format(1 - CalibSets.ALPHA, prec=0),
      calibrated_model.get_mean_survival_time_proj_interval(alpha=CalibSets.ALPHA))
