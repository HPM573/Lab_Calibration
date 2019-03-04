import CalibrationClasses as Cls
import CalibrationSettings as P
import SimPy.FigureSupport as Fig


# initialize a calibrated model
calibrated_model = Cls.CalibratedModel(csv_file_name='CalibrationResults.csv')
# simulate the calibrated model
calibrated_model.simulate(num_of_simulated_cohorts=P.NUM_SIM_COHORTS,
                          cohort_size=P.SIM_POP_SIZE,
                          time_steps=P.TIME_STEPS)

# plot the histogram of mean survival time
Fig.graph_histogram(
    data=calibrated_model.multiCohorts.multiCohortOutcomes.meanSurvivalTimes,
    title='Histogram of Mean Survival Time',
    x_label='Mean Survival Time (Year)',
    y_label='Count',
    bin_width=0.25,
    x_range=[7, 13])

# report mean and projection interval
print('Mean survival time and {:.{prec}%} projection interval:'.format(1 - P.ALPHA, prec=0),
      calibrated_model.get_mean_survival_time_proj_interval(alpha=P.ALPHA))
