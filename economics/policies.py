# This file defines the policies that can be modeled by the economic / outcomes model.
# It will define inputs used in each policy case used too calculate impact

# Policies are defined by a dict. Examples:
Policy_1 = dict(Economy_Open_Case_Number=1000, Hire_Interval=90, Trace=True, Test=False)
Policy_2 = dict(Economy_Open_Case_Number=1000, Hire_Interval=90, Trace=True, Test=True, Pop_test_pct='Adaptive')
Policy_3 = dict(Economy_Open_Case_Number=1000, Hire_Interval=90, Trace=True, Test=True, Pop_test_pct=1)

# Across Cases: Deaths, direct cost of hospitalization/ICU, economic cost of sickness, hospitalizations and deaths.
# Quarantine costs?

# Total deaths are output directly from the model. Peak Number ill is also output directly.


# Policy - Open Without Mitigation: No interventions, no (additional) distancing.

# Economy: Significant transmission will suppress economic output while ongoing.


# Policy - Open when diagnosed cases under X, with full contact tracing.

# Economy: Closed economy for time while cases > X, resumed economy when over X.

# Costs: Tracing costs. Assume workers are contracted for 3 months.
# Startup costs are fixed.
# In each 3-month period, we hire max(number of traces during period) * safety margin (10%)

# Policy - Open when cases under X, with full testing and contact tracing.

# Economy: Closed  economy for time while cases > X, resumed economy when over X.

# Costs: Tracing and testing costs. Assume workers are contracted for 3 months.
# In each 3-month period, we hire max(number of traces during period) * safety margin (10%)
# We do full population testing
# In each period, we also have labs + workers sufficient for max(number of tests during period) * safety margin (10%)