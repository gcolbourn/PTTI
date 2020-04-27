# This file lists the inputs for the economic cost model.

# It is still preliminary, and further inputs need to be moved from the spreadsheet.

# Parameters
UK_Population = 67886011
UK_GDP_Monthly = 186000000000

# NHS Costs
NHS cost of death
n_d,0

NHS cost of ICU case
n_i,0
£29,335
NHS cost of hospital (non-ICU) case
n_w,0
£2,422
Productivity cost of lost day of work
p_0
£119
Productivity cost of death
p_d,0
£358
Productivity cost of ICU case
p_i,0
£2,509
Productivity cost of hospital (non-ICU) case
p_w,0
£1,195
Productivity cost of non-hospital case
p_h,0
£358



# Fixed_costs


# Variable Tracing Costs
Number_of_Tracers =  UK_Population/1000
Cost_per_Tracer = 80 # Daily
Number_of_Tracing_Supervisors = Number_of_Tracers/50
Cost_per_Supervisor = 160 # Daily
Number_of_Tracing_Team_Leads = 343
Cost_per_Team_Lead = 300 # Daily

Phone_Credits =  Number_of_Tracers + Number_of_Tracing_Supervisors + Number_of_Tracing_Team_Leads
Phone_Credit_Costs = 5 # Daily

Smart_Phones = Number_of_Tracers/10 # For any tracers without phones. Minimal.

Rural_Pct = 0.17
Daily_Travel_Cost = 10
Travelers = (Number_of_Tracers * Rural_Pct) + Number_of_Tracing_Supervisors + Number_of_Tracing_Team_Leads

# Testing Costs

# PCR Machine Capabilities

PCR_Machines_Per_Lab = 10
Shifts_per_Day = 2
Hours_per_Shift = 9
Time_per_Batch = 0.5 # Half Hour batches
Tests_per_Batch = 96
Tests_per_Machine_per_Hour =  Tests_per_Batch / Time_per_Batch
Tests_per_Machine_per_Shift = Tests_per_Machine_per_Hour * Hours_per_Shift
Tests_per_Machine_per_Day = Tests_per_Machine_per_Shift * Shifts_per_Day

Lab_Techs_Per_Machine_Per_Shift = 2 # One to run the test, one to fill the wells.

Maximum_Daily_Tests = 10000000 # Per assumption. This will change based on model runs.
Total_tests = 3120000000 # Max * 312 Days. <- Should be Replaced by Modeled value.

#Setup / Fixed Costs
Needed_PCR_Machines = Maximum_Daily_Tests / Tests_per_Machine_per_Day
Needed_Labs = Needed_PCR_Machines / PCR_Machines_Per_Lab

#Personnel Costs
Lab_supervisors	= Needed_Labs
Max_Lab_Techs = Needed_PCR_Machines * Lab_Techs_Per_Machine_Per_Shift * Shifts_per_Day

Lab_staff_trainings = Max_Lab_Techs + Lab_supervisors # Retrainings every 3 months?


Lab_Tech_Salary = 200 # Per Shift
Lab_Supervisor_Salary = 300 # Per Day


# Variable Costs
Tests_Day_N = Maximum_Daily_Tests # This will be replaced with the daily number.

Lab_Techs_Day_N = Tests_Day_N / (Tests_per_Machine_per_Shift * Lab_Techs_Per_Machine_Per_Shift)
Labs_for_Day_N =  Tests_Day_N / (Tests_per_Machine_per_Day * PCR_Machines_Per_Lab)

Cost_per_Test_Kit = 4 # 3.50 for the testing, 0.50 for the home test kit.

Cost_Day_N =

