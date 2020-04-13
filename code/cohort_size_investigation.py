from single_cohort import *
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

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

##POPULATION ARRANGEMENT##
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


c_P = 0.01# relative force of infection from general population I1 and IA to nCovC patients.
c_H = 0.01# relative force of infection from general population infected to HCW.
c_PP = 0.5# 
c_HP = 2
c_PH = 2 #relative force of infection from patients to HCWs in cohort compared to in general public
c_HH = 0.5 # relative force of infection between HCWs in same cohort


test_incoming_factor = 0.95 #the probability an exposed/infected patient will be identified and not admitted.
       #if 1, then no infected individuals get through.  If 0, then there is no filter (except for symptomatic infections).

omega = 0.05
gamma_Q = 1./14  #departure rate from quarantine
ave_nonCOVIDduration = 14  #average hospital stay for non-COVID patient
#print('fix line above')

T=240 #end time


r''' We're ready to do the calculation
'''

H = {}
for NumberPatients in [50, 100, 200, 400, 800, 1600]:
    print(NumberPatients)
    if NumberPatients not in H:
        H[NumberPatients] =[]
    while len(H[NumberPatients])<50:
        if len(H[NumberPatients])%10==0:
            print('.')
        ts, pop_state, pop_FOI, Health_Facility_States, PFOIs, HFOIs = simulate(
            lambda_1, lambda_2, lambda_A, omega, gamma_E, gamma_I1, gamma_I2, 
            gamma_IA, gamma_Q, q, c_H, c_P,c_PP, c_HP, c_PH, c_HH, NumberPatients, 
            T, max_dt, ave_nonCOVIDduration, PatientsPerHCW=4, test_incoming_factor=test_incoming_factor)
        myCohort_state = np.array(Health_Facility_States)
        H[NumberPatients].append(myCohort_state[:,10][-1])
           
    fig, ax = plt.subplots(figsize=(8,4))
    n_bins = NumberPatients//4
    ax.hist(H[NumberPatients], n_bins, density=True, histtype = 'step', cumulative=True, label = NumberPatients)
    ax.axis(xmin=0)
    plt.savefig('{}_total_HCWs.png'.format(NumberPatients))
              
                                      
