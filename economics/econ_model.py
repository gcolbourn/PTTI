from economics.inputs import *
from economics.policies import *
from math import ceil, exp

# Right now, output['I'] is cumulative infections. TODO: Change this! (Pre-transform to that, maybe? Or fix logic?)

def test_output(max, growthrate, max_growth_day, days):
    model_output = dict()
    model_output['Days'] = days
    model_output['I'] = [max/(1+exp(-1*growthrate*(x-max_growth_day+10)) ) for x in range(0, days)]
    model_output['Total Deaths'] = model_output['I'][days-1]*0.01
    return(model_output)

test_output_0 = test_output(900, 0.25, 30, 110) # Good enough for a basic test.

from economics.policies import Test_Policy_1

def econ_outputs(policy, model_output):
    Days = model_output['Days'] # Total number of days modeled
    Costs = 0

    Trace_Outputs = dict()
    if policy['Trace']:
        Sick_Next = model_output['I'].copy()
        Sick_Next.insert(0, 0)
        Sick_Next = Sick_Next[0:-1] #This is offset by one to take the difference.
        # NOTE: Model output in compartment I is total sick, new_sick is daily difference, i.e. new people to trace.
        New_Sick = [a - b for a, b in zip(model_output['I'], Sick_Next)]


        if not policy['Test']:
            New_Sick = [num * Pct_Symptomatic for num in New_Sick]  # Only Sick people, since others aren't known to trace. (Now 50%)
            # This assumes we find all symptomatic people. TODO - We should fix that.

        # Now we need to find the number of tracers needed in each window.
        Tracers_Needed_Per_Hiring_Window = list()
        Period_Lengths = list()
        Windows = ceil(Days/policy['Hire_Interval'])
        for i in range(0, Windows):
            Period_Start = i*policy['Hire_Interval']
            Period_End = min((i+1)*policy['Hire_Interval'], Days)
            Period_Lengths.append(Period_End-Period_Start) #The last period might be shorter.
            Tracers_Needed_Per_Hiring_Window.append(max(New_Sick[Period_Start: Period_End])*Tracers_Per_Infected_Person)
            # We assume that each day, all new infected people have all contacts traced by a team.
            # This means to keep tracing on a one-day lag, i.e. by end of the next day,  we need enough people to trace
            # the maximum expected number in any one day during that period.

        #elif policy['Pop_test_pct'] <= 1:
        #else:
        #    raise NotImplementedError

        # The above now does a bunch of checking to get the numbers for tracers needed. Now we take the sum.
        Tracers_fixed_cost = Number_of_Tracing_Supervisors * Cost_per_Supervisor * Days + \
            Number_of_Tracing_Team_Leads * Cost_per_Team_Lead * Days + \
            max(Tracers_Needed_Per_Hiring_Window)

        Costs += Tracers_fixed_cost

        Tracers_cost = sum([(Cost_per_Tracer * a * b) for a, b in zip(Tracers_Needed_Per_Hiring_Window, Period_Lengths)])

        Costs += Tracers_cost

        Trace_Outputs=dict()
        Trace_Outputs['Tracers_Needed_Per_Hiring_Window'] = Tracers_Needed_Per_Hiring_Window
        Trace_Outputs['Period_Lengths'] = Period_Lengths
        Trace_Outputs['Tracers_fixed_cost'] = Tracers_fixed_cost
        Trace_Outputs['Tracers_cost'] = Tracers_cost

    if policy['Test']:
        if policy['Pop_test_pct'] == 'Adaptive':
            Daily_tests = model_output['Tests'] # Check test needs over time. Model gives me tests performed each day.
            period_length = policy['Hire_Interval'] #We need to split into periods of length

    return list(model_output['Total Deaths'], Costs, Trace_Outputs)
