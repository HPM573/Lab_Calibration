import numpy as np

import calib_sets as sets
import deampy.plots.histogram as hist
import deampy.plots.sample_paths as path
from MultiSurvivalModelClasses import MultiCohort

# find values of mortality probability at which the posterior should be evaluated
np.random.seed(1)
mortality_samples = np.random.uniform(
    low=sets.PRIOR_RANGE['Mortality Probability'][0],
    high=sets.PRIOR_RANGE['Mortality Probability'][1],
    size=sets.N_SAMPLES)

# create multiple cohorts
multiCohort = MultiCohort(
    ids=range(sets.N_SAMPLES),
    pop_sizes=[sets.SIM_POP_SIZE] * sets.N_SAMPLES,
    mortality_probs=mortality_samples  # [p1, p2, ....]
)

# simulate all cohorts
multiCohort.simulate(
    n_time_steps= sets.TIME_STEPS,
    seeds=range(sets.N_SAMPLES))

# plot the sample paths
path.plot_sample_paths(
    sample_paths=multiCohort.multiCohortOutcomes.survivalCurves,
    title='Survival Curves',
    x_label='Time-Step (Year)',
    y_label='Number Survived',
    transparency=0.25,
    x_range=(0, 100),
    file_name='figs/uncalibrated/survival_curves.png')

# plot the histogram of average survival time
hist.plot_histogram(
    data=multiCohort.multiCohortOutcomes.meanSurvivalTimes,
    title='Histogram of Mean Survival Time',
    x_label='Mean Survival Time (Year)',
    y_label='Count',
    x_range=[0, 25],
    bin_width=0.5,
    file_name='figs/uncalibrated/survival_times.png')

# create the histogram of the mortality probabilities
hist.plot_histogram(
    data=mortality_samples,
    title='Histogram of Mortality Probabilities',
    x_label='Mortality Probability',
    y_label='Counts',
    x_range=sets.PRIOR_RANGE['Mortality Probability'],
    bin_width=0.005,
    file_name='figs/uncalibrated/mortality_probs.png')

# print projected mean survival time (years)
print('Projected mean survival time (years)',
      multiCohort.multiCohortOutcomes.statMeanSurvivalTime.get_mean())

# print projection interval
print('95% projection (prediction, percentile, or uncertainty) interval of average survival time (years)',
      multiCohort.multiCohortOutcomes.statMeanSurvivalTime.get_PI(alpha=0.05))
