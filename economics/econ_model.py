from economics.inputs import *
from economics.policies import *

def econ_outputs(policy, model_output):
    Days = model_output['Days'] # Total number of days modeled
    Costs = 0

    Trace_Outputs = dict()
    if policy['Trace']:
        Tracers = [Max_Number_of_Tracers for i in range(0, Days)]  #Number Defined for each day
        Trace_Outputs['Daily_Tracers'] = Tracers
        # The above needs to depend on a bunch of checking.
        Tracers_fixed_cost = Number_of_Tracing_Supervisors * Cost_per_Supervisor * Days + \
            Number_of_Tracing_Team_Leads * Cost_per_Team_Lead * Days
        Tracer_cost = sum(Tracers * Cost_per_Tracer) + Number_of_Tracing_Supervisors * Cost_per_Supervisor + \
        Number_of_Tracing_Team_Leads * Cost_per_Team_Lead
        Costs = Costs + Tracer_cost

    if policy['Test']:
        if policy['Pop_test_pct'] == 'Adaptive':
            Daily_tests = model_output['Tests'] # Check test needs over time. Model gives me tests performed each day.
            period_length = policy['Hire_Interval'] #We need to split into periods of length


    return list(model_output['Total Deaths'],Costs, Trace_Outputs)
