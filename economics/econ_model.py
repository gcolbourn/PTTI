from economics.inputs import *
# from economics.policies import *
from math import ceil, exp
import yaml

Test = True
# Run model for testing:
if Test:
    from SEIRT.seirct_ode import SEIRCTODEMem
    Test_YAML = r'economics\ptti-discussion.yaml'
    with open(Test_YAML) as file:
        scenario = yaml.load(file, Loader=yaml.FullLoader)

    Model = SEIRCTODEMem()
    model_outputs = Model.run(t0=0, tsteps=365, tmax=365, state=Model.initial_conditions(**scenario['initial']))


def Econ_Outputs(model_outputs, YAML_file):
    with open(YAML_file) as file:
        scenario = yaml.load(file, Loader=yaml.FullLoader)

    Output = dict()
    Output['policy'] = scenario
    Output['model_output'] = model_outputs
    Output['timesteps'] = Output['model_output'][0]
    Days = len(model_outputs[0]) # Total number of days modeled
    Output['compartments'] = model_outputs[1]

    Costs = 0

    Output['econ'] = dict()
    Output['parameters'] = dict()
    Output['variables'] = dict()
    # GEt parameter values, which can change over time due to interventions.
    for var in scenario['parameters'].keys():
        Output['variables'][var] = [scenario['parameters'][var] for d in range(Days)] # Initial value.
    for intervention in scenario['interventions']:
        for var in intervention['parameters']:
            Output['variables'][var][intervention['time']:] = [intervention['parameters'][var] for i in range(Days-int(intervention['time']))]

    Output['I'] = [a + b for a, b in zip(Output['compartments'][:, 4], Output['compartments'][:,5])]


    Output['econ']['trace'] = [0 for d in range(Days)] # Number of people that must be traced.
    for d in range(Days):
        if Output['variables']['theta'][d] > 0:
            Output['econ']['trace'][d] = Output['variables']['c'][d] * Output['variables']['theta'][d] / (
                Output['variables']['gamma'][d] + Output['variables']['theta'][d] *
                (1 + Output['variables']['eta'][d] * Output['variables']['chi'][d]) ) * Output['compartments'][d, 6] #IU
    # if True: #policy['Trace']: # Currently, we always assumes we pay to trace,

    To_Trace = Output['econ']['trace']

    # Now we need to find the number of tracers needed in each window.
    Tracers_Needed_Per_Hiring_Window = []
    Period_Lengths = []
    Period_Starts = []
    Period_Ends = []
    Period_Start = 0
    for p in Output['policy']['interventions']:
        Period_End = p['time']
        Period_Starts.append(Period_Start)
        Period_Ends.append(Period_End)
        Period_Lengths.append(Period_End - Period_Start)
        Tracers_This_Period = max(To_Trace[Period_Start:Period_End]) * Time_to_Trace
        Tracers_Needed_Per_Hiring_Window.append(Tracers_This_Period)
        Period_Start = Period_End # For next period.
    # Last Period:
    Period_Starts.append(Period_Start)
    Period_Ends.append(Days)
    Period_Lengths.append(Days - Period_Start)
    Tracers_Needed_Per_Hiring_Window.append(max(To_Trace[Period_Start:Days]) * Time_to_Trace)

        # We assume that each day, all people whos contacts can be traced have all contacts traced by a team.
        # This means to keep tracing on a one-day lag, i.e. by end of the next day,  we need enough people to trace
        # the maximum expected number in any one day during that period.
    Max_Tracers = max(Tracers_Needed_Per_Hiring_Window)

    # The above now does a bunch of checking to get the numbers for tracers needed. Now we take the sum.
    Total_Max_Tracing_Workers = Number_of_Tracing_Supervisors + Number_of_Tracing_Team_Leads + Max_Tracers
    Recruitment_Costs = Hiring_Cost * Total_Max_Tracing_Workers

    Tracers_fixed_cost = (Number_of_Tracing_Supervisors * Cost_per_Supervisor * Days) + \
                         (Number_of_Tracing_Team_Leads * Cost_per_Team_Lead * Days) + \
                         Recruitment_Costs + Tracer_Training_Course_Cost + \
                         Cost_Per_Extra_Phones_for_Workers * Total_Max_Tracing_Workers

    Costs += Tracers_fixed_cost

    Supervisor_Travel_costs = Daily_Travel_Cost * (Number_of_Tracing_Supervisors + Number_of_Tracing_Team_Leads) * \
        sum(Period_Lengths)

    Tracer_Days_Needed = sum([(a * b) for a, b in zip(Tracers_Needed_Per_Hiring_Window, Period_Lengths)])
    Tracing_Worker_Travel_Costs = Daily_Travel_Cost * Rural_Pct  # Rural Tracers Travel Costs
    Tracers_cost = (Cost_per_Tracer + Tracing_Worker_Travel_Costs) * Tracer_Days_Needed + Supervisor_Travel_costs

    Costs += Tracers_cost

    Costs += Tracing_Daily_Public_Communications_Costs * Days

    # That is all costs from the original cost spreadsheet for Tracing.

    Trace_Outputs=dict()
    Trace_Outputs['Tracers_Needed_Per_Hiring_Window'] = Tracers_Needed_Per_Hiring_Window
    Trace_Outputs['Period_Lengths'] = Period_Lengths
    Trace_Outputs['Tracers_fixed_cost'] = Tracers_fixed_cost
    Trace_Outputs['Tracers_cost'] = Tracers_cost

    #Testing Level is now set by Theta, as such:.
    # theta: 100000 / N
    Test_Outputs = dict()

    Testing_Costs = 0 # We will add to this.

    Output['econ']['tests'] = [0 for d in range(Days)] # Number of people tested.
    for d in range(Days):
        if Output['variables']['theta'][d] > 0:
            Output['econ']['tests'][d] = Output['policy']['initial']['N'] * Output['variables']['theta'][d]

    Daily_tests = Output['econ']['tests']  # Test needs over time.

    # Model gives me tests performed each day - or we need to calculate it here. TODO: Check this
    Max_Tests_In_Hiring_Window = []  # This determines how many lab techs we need per period.
    Labs_Per_Window = []
    Maximum_tests = max(Daily_tests)

    for i in range(0, len(Period_Lengths)):
        Window_Tests_max = max(Daily_tests[Period_Starts[i]: Period_Ends[i]])
        Max_Tests_In_Hiring_Window.append(Window_Tests_max)
        Machines_Needed_In_Window = Window_Tests_max / Tests_per_Machine_per_Day
        Labs_Per_Window.append(Machines_Needed_In_Window / PCR_Machines_Per_Lab) # Others can be shut down.

    # Fixed Testing Costs:
    Total_Required_PCR_Machines = max(Max_Tests_In_Hiring_Window) / (Tests_per_Machine_per_Day)
    Max_Lab_Staff = Total_Required_PCR_Machines / PCR_Machines_Per_Lab  # 1 Supervisor per lab.
    Max_Lab_Staff += Total_Required_PCR_Machines / (Shifts_per_Day * Lab_Techs_Per_Machine_Per_Shift)
    Testing_Costs += Total_Required_PCR_Machines*PCR_Machines_Cost # Buy Machines
    Testing_Costs += Max_Lab_Staff * Hiring_Cost # Hire Staff. One Time Cost.

    Total_Required_Tests = sum(Daily_tests)
    Testing_Costs += Total_Required_Tests * Cost_Per_PCR_Test

    for p in range(0, len(Period_Lengths)):
        Machines_In_Window = Max_Tests_In_Hiring_Window[p] / Tests_per_Machine_per_Day
        Labs_In_Window = Machines_In_Window / PCR_Machines_Per_Lab
        Curr_Period_Days = Period_Lengths[p]
        Testing_Costs += Lab_Overhead_Cost_Daily * Curr_Period_Days
        # Labs in use are all fully staffed.
        # Staff Costs
        Supervisors = Labs_In_Window  # 1 Supervisor per lab.
        Testing_Costs += Supervisors * Lab_Supervisor_Salary * Curr_Period_Days
        Daily_Lab_Workers = Labs_In_Window * PCR_Machines_Per_Lab * Shifts_per_Day * Lab_Techs_Per_Machine_Per_Shift
        Testing_Costs += Daily_Lab_Workers * Lab_Tech_Salary
        Testing_Costs += Staff_Training_Cost * (Daily_Lab_Workers + Supervisors) # Once Per Period

        Testing_Costs += Machines_In_Window * Curr_Period_Days * PCR_Machine_Daily_Maintenance

    Testing_Costs += sum(Daily_tests)*Cost_Per_PCR_Test # Cost for the actual tests.


    Test_Outputs = dict()
    Test_Outputs['Max_Labs'] = Total_Required_PCR_Machines / PCR_Machines_Per_Lab
    Test_Outputs['Max_Total_Staff'] = Max_Lab_Staff
    Test_Outputs['Total_Tests'] = sum(Daily_tests)

    Costs += Testing_Costs

    Output['Costs'] = Costs
    Output['Trace_Outputs'] = Trace_Outputs
    Output['Test_Outputs'] = Test_Outputs

    return Output  #  What else is useful to graph / etc?

#
# test_output_0 = test_output(900, 0.25, 30, 110)  # Good enough for a basic test.
#
# from economics.policies import Test_Policy_2
#
# econ_outputs(Test_Policy_2, test_output_0)