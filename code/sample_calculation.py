from single_cohort import *
import matplotlib.pyplot as plt
import numpy as np

#
###PARAMETERS###
max_dt = 0.001  #maximum time step

##TRANSMISSION RATES##
lambda_1 = 0.25    #transmission parameters for general population
lambda_2 = 2./7
lambda_A = lambda_1

##TRANSITION RATES##  (rates to transition out of each state)
gamma_E = 1./3   
gamma_I1 = 1./2
gamma_I2 = 1./7
gamma_IA = 1./9

q = 0.5  #probability of arriving into asymptomatic state rather than presymptomatic.

test_incoming_factor = 0 #the probability an infected patient will be identified and not admitted.

##POPULATION ARRANGEMENT##
NumberPatients = 1000 #initial number of nonCovC patients
#NumberCohorts = 1  #number of nCovC subcohorts
#    Nhat = NumberPatients/NumberCohorts  #typical subCohort size
#    b = Nhat/ave_nonCOVIDduration   #entry rate required to have the typical stay be ave_nonCOVIDduration and size be Nhat.
    
   
# Hx_count = int(round(NumberPatients/4))  # number of HCW in nonCov cohort Assuming 4 patients/HCW
#     
# myCohort = Cohort(Hx_count,NumberPatients)
# HCW_statuses = [Hx_count, 0, 0, 0, 0, 0]
# HCW_trans_rates = [] 
#     
# #P_statuses = Patient_nCovC_status = [np.array([int(NumberPatients/NumberCohorts), 0, 0]) for n in range(NumberCohorts)]
# Patient_trans_rates = []
# 
# Pop_status = np.array([0.99, 0.01, 0, 0, 0, 0])  #initial condition for population, normalized to 1]


##TRANSMISSION PARAMETERS##


c_P = 0.1# relative force of infection from general population I1 and IA to nCovC patients.
c_H = 0.1# relative force of infection from general population infected to HCW.
c_PP = 0.5# 
c_HP = 2
c_PH = 2 #relative force of infection from patients to HCWs in cohort compared to in general public
c_HH = 0.5 # relative force of infection between HCWs in same cohort

omega = 0*0.05
gamma_Q = 1./14  #departure rate from quarantine
ave_nonCOVIDduration = 14  #average hospital stay for non-COVID patient
#print('fix line above')

T=100 #end time


r''' We're ready to do the calculation
'''

ts, pop_state, pop_FOI, Health_Facility_States, PFOIs, HFOIs = simulate(
            lambda_1, lambda_2, lambda_A, omega, gamma_E, gamma_I1, gamma_I2, 
            gamma_IA, gamma_Q, q, c_H, c_P,c_PP, c_HP, c_PH, c_HH, NumberPatients, 
            T, max_dt, ave_nonCOVIDduration, PatientsPerHCW=4, test_incoming_factor=test_incoming_factor)
  
           
                    
                             
                                      
                                               
                                                        
###PRODUCING OUTPUTS                                                                 
                                                                                   
#Patient_transmission_rates = [sum(L) for L in Patient_nCovC_trans_rates]
#HCW_transmission_rates = [sum(L) for L in HCW_nCovC_trans_rates]
plt.figure(3)
plt.clf()
myCohort_state = np.array(Health_Facility_States)# for Cohortstate in Health_Facility_States]
pop_state = np.array(pop_state)
#plt.plot(ts, pop_state[:,1], label = 'population E')

plt.plot(ts, pop_state[:,0], label = 'susceptible general population')
plt.plot(ts, sum(pop_state[:,i] for i in [1,2,3,4]), label = 'infected general population')
#Cohort_state =sum(Cohortstates)

'''
The structure of pop_state is such that 
pop_state[:,0] is an array of 'S'
pop_state[:,1] is an array of 'E'
pop_state[:,2] is an array of I1
[:,3] -> I2
[:,4] -> IA
[:,5] -> R
'''


r'''  Cohort_state gives the states of the patients first and then the states 
of the HCWs.
By looking at the labels in the plotting commands you can figure out which is 
which.
'''

plt.plot(ts, myCohort_state[:,5]/myCohort_state[0,5], '-.', label = 'Susceptible  HCW' )
# plt.plot(ts, myCohort_state[:,6]/myCohort_state[0,5], '-.', label = 'E HCW')
# plt.plot(ts, myCohort_state[:,7]/myCohort_state[0,5], '-.', label = 'I1 HCW')
# plt.plot(ts, myCohort_state[:,8]/myCohort_state[0,5], '-.', label = 'IA HCW')
plt.plot(ts, sum(myCohort_state[:,j] for j in [6,7,8])/myCohort_state[0,5], '-.', label = 'Infected active HCW')
plt.plot(ts, myCohort_state[:,9]/myCohort_state[0,5], '-.', label = 'Quarantined HCW')
plt.plot(ts, myCohort_state[:,10]/myCohort_state[0,5], '-.', label = 'Recovered HCW')
plt.plot(ts, myCohort_state[:,0]/myCohort_state[0,0], '--', label = 'Susceptible Patient')
#plt.plot(ts, myCohort_state[:,1]/myCohort_state[0,0], '--', label = 'E  Patient')
#plt.plot(ts, myCohort_state[:,2]/myCohort_state[0,0], '--', label = 'I1  Patient')
#plt.plot(ts, myCohort_state[:,3]/myCohort_state[0,0], '--', label = 'IA  Patient')
plt.plot(ts, sum(myCohort_state[:,j] for j in [1,2,3])/myCohort_state[0,0], '-.', label = 'Infected Patient')
#plt.plot(ts, myCohort_state[:,4]/myCohort_state[0,0], '--', label = 'R  Patient')
#plt.plot(ts, (myCohort_state[:,0]+myCohort_state[:,1]+myCohort_state[:,2]+myCohort_state[:,3]+myCohort_state[:,4])/myCohort_state[0,0])
plt.legend()
plt.show()  

plt.figure(4)
plt.clf()
plt.plot(ts, pop_FOI, label = 'public FOI')
plt.plot(ts, PFOIs, '--', label = 'Patient FOI')
plt.plot(ts, HFOIs, '-.', label = 'HCW FOI')
plt.legend()
plt.show()

#print(HCW_nCovC_status)


plt.figure(3)
plt.savefig('base_case.png')
plt.figure(4)
plt.savefig('base_FOI.png')
            
