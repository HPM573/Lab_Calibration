import deampy.plots.histogram as hist
import deampy.plots.sample_paths as path
import numpy as np

import CalibrationSettings as Sets
import MultiSurvivalModelClasses as Cls

# find values of mortality probability at which the posterior should be evaluated
np.random.seed(1)
mortality_samples = np.random.uniform(
    low=Sets.PRIOR_L,
    high=Sets.PRIOR_U,
    size=Sets.PRIOR_N)

# create multiple cohorts
multiCohort = Cls.MultiCohort(
    ids=range(Sets.PRIOR_N),
    pop_sizes=[Sets.SIM_POP_SIZE] * Sets.PRIOR_N,
    mortality_probs=mortality_samples  # [p1, p2, ....]
)

# simulate all cohorts
multiCohort.simulate(Sets.TIME_STEPS)

# plot the sample paths
path.plot_sample_paths(
    sample_paths=multiCohort.multiCohortOutcomes.survivalCurves,
    title='Survival Curves',
    x_label='Time-Step (Year)',
    y_label='Number Survived',
    transparency=0.5)

# plot the histogram of average survival time
hist.plot_histogram(
    data=multiCohort.multiCohortOutcomes.meanSurvivalTimes,
    title='Histogram of Mean Survival Time',
    x_label='Mean Survival Time (Year)',
    y_label='Count',
    x_range=[2.5, 21.5],
    bin_width=2)

# create the histogram of the mortality probabilities
hist.plot_histogram(
    data=mortality_samples,
    title='Histogram of Mortality Probabilities',
    x_label='Mortality Probability',
    y_label='Counts',
    x_range=[Sets.PRIOR_L, Sets.PRIOR_U])

# print projected mean survival time (years)
print('Projected mean survival time (years)',
      multiCohort.multiCohortOutcomes.statMeanSurvivalTime.get_mean())

# print projection interval
print('95% projection (prediction, percentile, or uncertainty) interval of average survival time (years)',
      multiCohort.multiCohortOutcomes.statMeanSurvivalTime.get_PI(alpha=0.05))
