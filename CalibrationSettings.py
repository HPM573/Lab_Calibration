import scipy.stats as stat

# details of a clinical study estimating the mean survival time
OBS_N = 150        # number of patients involved in the study
OBS_MEAN = 10.2    # estimated mean survival time
OBS_HL = 1.5      # half-length
OBS_ALPHA = 0.05   # significance level

SIM_POP_SIZE = 500       # population size of simulated cohorts
TIME_STEPS = 1000        # length of simulation
ALPHA = 0.05             # significance level for calculating confidence intervals
NUM_SIM_COHORTS = 100   # number of simulated cohorts used to calculate prediction intervals

# the standard deviation of the mean survival time reported in the clinical study
# assumes that the reported confidence interval in this study is a t-confidence interval
OBS_STDEV = OBS_HL / stat.t.ppf(1 - OBS_ALPHA / 2, OBS_N-1)

# how to sample the posterior distribution of mortality probability
# minimum, maximum and the number of samples for the mortality probability
PRIOR_L, PRIOR_U, PRIOR_N = 0.05, 0.15, 100
