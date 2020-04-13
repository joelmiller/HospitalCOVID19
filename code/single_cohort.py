import matplotlib.pyplot as plt
import random
import numpy as np


class Cohort():
    '''This class handles a subcohort within the non-COVID-19 cohort.  
    It tracks how many people of each type there are.
     
    '''
    def __init__(self, NH, NP):
        self.SP = NP
        self.SH = NH
        self.NP = NP
        self.NH = NH
        
        self.EP = 0
        self.EH = 0
        
        self.I1P = 0
        self.I1H = 0
        
        self.I2P = 0
        self.I2H = 0
        
        self.IAP = 0
        self.IAH = 0
        
        self.QH = 0
        
        self.RP = 0
        self.RH = 0
        
        self.to_EP_rate = 0
        self.to_EH_rate = 0
        
        self.to_I1P_rate =0
        self.to_I1H_rate = 0
        
        self.to_I2P_rate = 0
        self.to_I2H_rate = 0
        
        self.to_IAP_rate = 0
        self.to_IAH_rate = 0
        
        self.to_RP_rate = 0
        self.to_RH_rate = 0
        
        
        self.departure_rate = 0
        self.arrival_rate = 0
        self.total_rate = 0
        
        self.P1taus = []
        self.PAtaus = []
        self.Ptau_sum = 0
        
        self.H1taus = []
        self.HAtaus = []
        self.Htau_sum = 0
        
        self.cumulative_rate = []
        
    def calculate_rates(self,  popS, popE, popI1, popIA, popR, poplambda, #myCohort, #I1oH, IAoH, 
                                CH, CP, CPP, CHxP, CPHx, CHxHx, b, Nhat, lambda1, 
                                lambdaA, omega, gammaE, gammaI1, gammaI2, 
                                gammaIA, gammaQ, q, test_incoming_factor=0):
                                
        assert(lambdaA == lambda1)
        self.total_rate = 0
        
        
        self.to_EP_rate = 0
        self.to_EH_rate = 0
        # N_H = sum(cohort.NH for cohort in cohorts if cohort != self)
        # if N_H>0:
        #     self.to_EH_rate = CHH* sum(lambda1*cohort.I1H +lambdaA*cohort.IAH for 
        #                            cohort in cohorts if cohort !=self)*self.SH/ N_H
        self.to_EH_rate += CH*poplambda  * self.SH
        # if CoHHx>0:
        #     CoHHx*(lambda1*I1oH+lambdaA*IAoH)/(N_H+self.NH)) * self.SH
        #self.to_EH_rate += (CH*(lambda1*popI1 + lambdaA*popIA)) *self.SH
        self.to_EP_rate += (CP*(lambda1*popI1 + lambdaA*popIA)) *self.SP
        
        #Htau_sum = sum(self.Htaus)
        #Ptau_sum = sum(self.Ptaus)
        #print(self.Htau_sum, sum(self.H1taus)+sum(self.HAtaus))
        self.to_EH_rate += (CPHx*self.Ptau_sum + CHxHx*self.Htau_sum)*self.SH/self.NH
        self.to_EP_rate += (CPP*self.Ptau_sum + CHxP*self.Htau_sum)*self.SP/self.NP
        
        self.total_rate += self.to_EH_rate + self.to_EP_rate
                
        self.to_I1P_rate = (1-q)*gammaE*self.EP 
        self.to_I2P_rate = gammaI1 * self.I1P
        self.to_IAP_rate = q*gammaE*self.EP
        self.from_IAP_rate = gammaIA*self.IAP
        self.removal_by_test_P_rate = omega * (self.EP + self.I1P + self.IAP)
        self.total_rate += self.to_I1P_rate + self.to_I2P_rate + self.to_IAP_rate + self.to_RP_rate + self.removal_by_test_P_rate
        
        self.to_I1H_rate = (1-q)*gammaE*self.EH
        #self.to_I2H_rate = gammaI1 * self.I1H
        #XXXX
        self.to_IAH_rate = q*gammaE*self.EH
        self.from_IAH_rate = gammaIA*self.IAH
        
        self.total_rate += self.to_I1H_rate + self.to_IAH_rate + self.from_IAH_rate

        self.to_Q_by_test_rate = omega * (self.EH + self.I1H + self.IAH )
        self.to_Q_by_symptom_rate = gammaI1*self.I1H
        self.from_Q_rate = gammaQ*self.QH

        self.total_rate += self.to_Q_by_test_rate + self.to_Q_by_symptom_rate + self.from_Q_rate
        
        arrival_weight = (popS+ popR)
        arrival_weight += (1-test_incoming_factor)*(popE+popI1+popIA)
            
        self.Sarrival_rate = b*popS/arrival_weight
        self.Rarrival_rate = b*popR/arrival_weight
        
        
        self.Earrival_rate = (1-test_incoming_factor)*b*popE/arrival_weight
        self.I1arrival_rate = (1-test_incoming_factor)*b*popI1/arrival_weight
        self.IAarrival_rate = (1-test_incoming_factor)*b*popIA/arrival_weight
        self.total_rate += self.IAarrival_rate + self.Earrival_rate + self.I1arrival_rate
        self.total_rate += self.Sarrival_rate + self.Rarrival_rate 

        self.Sdeparture_rate = b*(self.SP)/Nhat
        self.Edeparture_rate = b*self.EP/Nhat
        self.I1departure_rate = b*self.I1P/Nhat
        self.IAdeparture_rate = b*self.IAP/Nhat
        self.Rdeparture_rate = b*self.RP/Nhat
        
        self.total_rate += self.Sdeparture_rate + self.Edeparture_rate 
        self.total_rate += self.I1departure_rate + self.IAdeparture_rate + self.Rdeparture_rate
        
        cr = []
        cr.append(self.to_EP_rate)
        cr.append(cr[-1] + self.to_I1P_rate)
        cr.append(cr[-1] + self.to_I2P_rate)
        cr.append(cr[-1] + self.to_IAP_rate)
        cr.append(cr[-1] + self.from_IAP_rate)
        cr.append(cr[-1] + self.removal_by_test_P_rate)
        
        cr.append(cr[-1] + self.to_EH_rate)
        cr.append(cr[-1] + self.to_I1H_rate)
        #cr.append(cr[-1] + self.to_I2H_rate)
        cr.append(cr[-1] + self.to_IAH_rate)
        #cr.append(cr[-1] + self.from_I2H_rate)
        cr.append(cr[-1] + self.from_IAH_rate)
        cr.append(cr[-1] + self.to_Q_by_test_rate)
        cr.append(cr[-1] + self.to_Q_by_symptom_rate)
        cr.append(cr[-1] + self.from_Q_rate)
        
        cr.append(cr[-1] + self.Sarrival_rate)
        cr.append(cr[-1] + self.Earrival_rate)
        cr.append(cr[-1] + self.I1arrival_rate)
        cr.append(cr[-1] + self.IAarrival_rate)
        cr.append(cr[-1] + self.Rarrival_rate)
        cr.append(cr[-1] + self.Sdeparture_rate)
        cr.append(cr[-1] + self.Edeparture_rate)
        cr.append(cr[-1] + self.I1departure_rate)
        cr.append(cr[-1] + self.IAdeparture_rate)
        cr.append(cr[-1] + self.Rdeparture_rate)
        
        #print(self.arrival_rate, self.Sdeparture_rate,self.SP, self.EP, self.I1P, self.IAP, self.RP, Nhat)
        #print(self.total_rate, cr[-1])
        self.cumulative_rate = cr
        #print(self.to_EH_rate,'EH')
        self.new_tauP = lambda : lambdaA
        
        self.new_tauH = lambda : lambdaA
        
    def do_action(self):
        #print(self.NH, self.EP)
        #print(self.cumulative_rate)
        r= random.random()*self.total_rate
        #print(r, self.total_rate, self.cumulative_rate)
        if r< self.cumulative_rate[0]:        #cr.append(self.to_EP_rate)
            self.SP -= 1
            self.EP += 1
        elif r< self.cumulative_rate[1]:   #cr.append(cr[-1] + self.to_I1P_rate)
            self.EP -= 1
            #print('I1P reduce EP')
            self.I1P += 1
            self.P1taus.append(self.new_tauP())
            self.Ptau_sum += self.P1taus[-1]
        elif r< self.cumulative_rate[2]:     #cr.append(cr[-1] + self.to_I2P_rate) 
            self.I1P -= 1
            self.NP -=1
            index = random.randrange(len(self.P1taus))
            self.P1taus[index], self.P1taus[-1] = self.P1taus[-1], self.P1taus[index] 
            self.Ptau_sum -= self.P1taus.pop()
        elif r<self.cumulative_rate[3]:        #cr.append(cr[-1] + self.to_IAP_rate)
            #print('IAP reduce EP')
            self.EP -= 1
            self.IAP += 1
            self.PAtaus.append(self.new_tauP())
            self.Ptau_sum += self.PAtaus[-1]
        elif r<self.cumulative_rate[4]:        #cr.append(cr[-1] + self.from_IAP_rate)
            self.IAP -= 1
            self.RP += 1
            index = random.randrange(len(self.PAtaus))
            self.PAtaus[index], self.PAtaus[-1] = self.PAtaus[-1], self.PAtaus[index] 
            self.Ptau_sum -= self.PAtaus.pop()     
        elif r<self.cumulative_rate[5]: #cr.append(cr[-1] + self.removal_by_test_P_rate)
            inf_P_size = self.EP + self.IAP + self.I1P
            r2 = random.random()*inf_P_size
            self.NP -= 1
            if r2 < self.EP:
                self.EP -=1
            elif r2 < self.EP + self.IAP:
                self.IAP -= 1
                index = random.randrange(len(self.PAtaus))
                self.PAtaus[index], self.PAtaus[-1] = self.PAtaus[-1], self.PAtaus[index] 
                self.Ptau_sum -= self.PAtaus.pop()     
            else:#if r2 < self.EP + self.IAP + self.I1P:
                self.I1P -= 1
                index = random.randrange(len(self.P1taus))
                self.P1taus[index], self.P1taus[-1] = self.P1taus[-1], self.P1taus[index] 
                self.Ptau_sum -= self.P1taus.pop()     
                
        elif r<self.cumulative_rate[6]:        #cr.append(cr[-1] + self.to_EH_rate)
            #print('H infection')
            self.SH -= 1
            self.EH += 1
        elif r< self.cumulative_rate[7]:        #cr.append(cr[-1] + self.to_I1H_rate)
            self.EH -=1
            self.I1H += 1
            self.H1taus.append(self.new_tauH())
            self.Htau_sum += self.H1taus[-1]
        # elif r< self.cumulative_rate[7]:        #cr.append(cr[-1] + self.to_I2H_rate)
        #     self.I1H -= 1
        #     self.I2H += 1
        #     self.NH -= 1
        #     index = random.randrange(len(self.H1taus))
        #     self.H1taus[index], self.H1taus[-1] = self.H1taus[-1], self.H1taus[index] 
        #     self.Htau_sum -= self.H1taus.pop()
        elif r< self.cumulative_rate[8]:        #cr.append(cr[-1] + self.to_IAH_rate)
            self.EH -= 1
            self.IAH += 1
            self.HAtaus.append(self.new_tauH())
            self.Htau_sum += self.HAtaus[-1]            
        # elif r<self.cumulative_rate[9]:        #cr.append(cr[-1] + self.from_I2H_rate)
        #     self.RH += 1
        #     self.I2H -= 1
        #     self.NH +=1
        elif r< self.cumulative_rate[9]:        #cr.append(cr[-1] + self.from_IAH_rate)
            self.RH += 1
            self.IAH -=1
            index = random.randrange(len(self.HAtaus))
            self.HAtaus[index], self.HAtaus[-1] = self.HAtaus[-1], self.HAtaus[index] 
            self.Htau_sum -= self.HAtaus.pop()            

        elif r< self.cumulative_rate[10]:        #cr.append(cr[-1] + self.to_Q_by_test_rate)
            inf_H_size = self.EH + self.IAH + self.I1H
            r2 = random.random()*inf_H_size
            self.QH += 1
            self.NH -= 1
            if r2 < self.EH:
                self.EH -=1
            elif r2 < self.EH + self.IAH:
                self.IAH -= 1
                index = random.randrange(len(self.HAtaus))
                self.HAtaus[index], self.HAtaus[-1] = self.HAtaus[-1], self.HAtaus[index] 
                self.Htau_sum -= self.HAtaus.pop()            
            else:#if r2 < self.EP + self.IAP + self.I1P:
                self.I1H -= 1
                index = random.randrange(len(self.H1taus))
                self.H1taus[index], self.H1taus[-1] = self.H1taus[-1], self.H1taus[index] 
                self.Htau_sum -= self.H1taus.pop()            
        elif r< self.cumulative_rate[11]:     #cr.append(cr[-1] + self.to_Q_by_symptom_rate)
            self.NH -=1
            self.QH += 1
            self.I1H -= 1
            index = random.randrange(len(self.H1taus))
            self.H1taus[index], self.H1taus[-1] = self.H1taus[-1], self.H1taus[index] 
            self.Htau_sum -= self.H1taus.pop()            
        elif r<self.cumulative_rate[12]:   #cr.append(cr[-1] + self.from_Q_rate)
            self.NH += 1
            self.RH += 1
            self.QH -= 1
            
        elif r< self.cumulative_rate[13]:        #cr.append(cr[-1] + self.Sarrival_rate)
            self.SP += 1
            self.NP += 1
        elif r< self.cumulative_rate[14]:        #cr.append(cr[-1] + self.Earrival_rate)
            self.EP += 1
            self.NP += 1
        elif r< self.cumulative_rate[15]:        #cr.append(cr[-1] + self.I1arrival_rate)
            self.I1P += 1
            self.NP += 1
            self.P1taus.append(self.new_tauP())
            self.Ptau_sum += self.P1taus[-1]
        elif r< self.cumulative_rate[16]:        #cr.append(cr[-1] + self.IAarrival_rate)
            self.IAP += 1
            self.NP += 1
            self.PAtaus.append(self.new_tauP())
            self.Ptau_sum += self.PAtaus[-1]
        elif r< self.cumulative_rate[17]:        #cr.append(cr[-1] + self.Rarrival_rate)
            self.RP += 1
            self.NP += 1
        elif r< self.cumulative_rate[18]:        #cr.append(cr[-1] + self.Sdeparture_rate)
            self.SP -= 1
            self.NP -=1
        elif r<self.cumulative_rate[19]:        #cr.append(cr[-1] + self.Edeparture_rate)
            self.EP -=1
            #print('Edepart reduce EP', self.Edeparture_rate)
            self.NP -=1
        elif r< self.cumulative_rate[20]:        #cr.append(cr[-1] + self.I1departure_rate)
            self.I1P -= 1
            self.NP -=1
            index = random.randrange(len(self.P1taus))
            self.P1taus[index], self.P1taus[-1] = self.P1taus[-1], self.P1taus[index] 
            self.Ptau_sum -= self.P1taus.pop()            
        elif r< self.cumulative_rate[21]:        #cr.append(cr[-1] + self.IAdeparture_rate)
            self.IAP -= 1
            self.NP -= 1
            #print(self.IAP)
            index = random.randrange(len(self.PAtaus))
            self.PAtaus[index], self.PAtaus[-1] = self.PAtaus[-1], self.PAtaus[index] 
            self.Ptau_sum -= self.PAtaus.pop()            
        elif r< self.cumulative_rate[22]:         #cr.append(cr[-1] + self.Rdeparture_rate)
            self.RP -=1
            self.NP -=1
        else:
            print(r, self.cumulative_rate)
            
        assert(self.I1H == len(self.H1taus))
        #print(self.IAH, len(self.HAtaus))
        assert(self.IAH == len(self.HAtaus))
        assert(self.I1P == len(self.P1taus))
        assert(self.IAP == len(self.PAtaus))
        assert(self.NH== self.SH + self.EH + self.IAH + self.I1H + self.RH)
        assert(self.NP ==self.SP + self.EP + self.IAP + self.I1P + self.RP)
        #print(self.Edeparture_rate, self.to_Q_by_test_rate)
        #print('XXX',len(self.cumulative_rate), self.cumulative_rate[-1], self.total_rate)
        
def dPublicdt(X, parameters):
    '''
    Used to give the derivative of the proportion of the public in each state.
    '''
    S, E, I1, I2, IA, R = X
    lambda1, lambda2, lambdaA, gammaE, gammaI1, gammaI2, gammaIA, q = parameters
    #print(X, '\n', parameters)
    #1/0
    FOI= (lambda1*I1 + lambda2*I2 + lambdaA*IA)
    dS = -FOI*S
    dE = -dS - gammaE*E
    dI1 = (1-q)*gammaE*E - gammaI1*I1
    dI2 = gammaI1*I1 - gammaI2*I2
    dIA = q*gammaE*E - gammaIA*IA
    dR = gammaI2 *I2 + gammaIA*IA
    return np.array([dS, dE, dI1, dI2, dIA, dR]), FOI

# def dCovC_HCW(X, parameters):
#     '''
#     This handles the derivatives for the HCWs who are treating COVID-19 patients
#     '''
#     S, E, I1, I2, IA, R = X
#     CH, lamb, K, CHH, lambda1, lambdaA, FoHH , gammaE, gammaI1, gammaI2, gammaIA, q= parameters
#     #Chlambda is CH * lambda, and FoHH =CoHH*lambda1*sum_j I1H_j is force of infection from other cohort 
#     FOI = 0#(CH*lamb + K + CHH * (lambda1*I1+lambdaA*IA)/(S+E+I1+IA+R) + FoHH/(S+E+I1+IA+R))
#     dS = -FOI*S
#     dE = -dS - gammaE * E
#     dI1 = (1-q)*gammaE * E - gammaI1*I1
#     dI2 = gammaI1*I1 - gammaI2*I2
#     dIA = q*gammaE *E - gammaIA*IA
#     dR = gammaI2*I2 + gammaIA*IA
#     #print(dS)
#     #print(CH*lamb, K, CHH, lambda1*I1/(S+E+I1), FoHH)
#     return 0*np.array([dS, dE, dI1, dI2, dIA, dR]), 0*FOI
#     



def simulate(lambda_1, lambda_2, lambda_A, omega, gamma_E, gamma_I1, gamma_I2, gamma_IA, gamma_Q,  q, 
            c_H, c_P, c_PP, c_HxP, c_PHx, c_HxHx, NumberPatients, T, 
            max_dt, ave_nonCOVIDduration, PatientsPerHCW, test_incoming_factor = 0):
            
    #oH_count = 0#initial number of CovC health care workers
    #oH_status = [oH_count, 0, 0, 0, 0, 0]  #S, E, I1, I2, IA, Rfor Covid-19 cohort of HCW
    K=0 
    #x_oHHx  
    
    Nhat = NumberPatients#/NumberCohorts  #typical subCohort size
    b = Nhat/ave_nonCOVIDduration   #entry rate required to have the typical stay be ave_nonCOVIDduration and size be Nhat.
    
    Hx_count = int(round(NumberPatients/PatientsPerHCW))  # number of HCW in nonCov cohort Assuming 4 patients/HCW
    myCohort = Cohort(Hx_count,NumberPatients)
    #HCW_statuses = [int(Hx_count/NumberCohorts), 0, 0, 0, 0, 0] 
    #HCW_nCovC_trans_rates = [[] for n in range(NumberCohorts)]
    
    #P_statuses = Patient_nCovC_status = [np.array([int(NumberPatients/NumberCohorts), 0, 0]) for n in range(NumberCohorts)]
    #Patient_nCovC_trans_rates = [[] for n in range(NumberCohorts)]
    
    Pop_status = np.array([0.99, 0.01, 0, 0, 0, 0])  #initial condition for population, normalized to 1]

    t = 0

    pop_state =[]
    pop_FOI = []
    #oH_state = []
    #oH_FOI = []
    Health_Facility_States = [] 
    HCW_FOIs = [] 
    P_FOIs = [] 
    ts = []
    while t<T:
        #print(t)
        S, E, I1, I2, IA, R = Pop_status
        #SoH, EoH, I1oH, I2oH, IAoH, RoH= oH_status
        pop_state.append(Pop_status)
        #oH_state.append(oH_status)
        #for nonCOVstate, sCohort in zip(Health_Facility_States, subCohorts):
        Health_Facility_States.append([myCohort.SP, myCohort.EP, myCohort.I1P, myCohort.IAP, myCohort.RP, myCohort.SH, myCohort.EH, myCohort.I1H, myCohort.IAH, myCohort.QH, myCohort.RH])
        ts.append(t)
        poplambda = lambda_1 * I1 + lambda_2*I2 + lambda_A*IA
        
        #for sCohort in subCohorts:
        myCohort.calculate_rates(S, E, I1, IA, R, poplambda, #myCohort, #I1oH, IAoH, 
                                    c_H, 
                                    c_P, c_PP, c_HxP, c_PHx, c_HxHx, b, Nhat, lambda_1, 
                                    lambda_A, omega, gamma_E, gamma_I1, gamma_I2, 
                                    gamma_IA, gamma_Q, q, test_incoming_factor)
            #print('NH', sCohort.NH)
        
        total_rate = myCohort.total_rate# for sCohort in subCohorts)
    
        #This is where we do a Gillespie-like step.  We assume the general population
        #(which follows deterministic dynamics) is unchanging until the next Gillespie 
        #step, but we don't allow large steps.  If that step is bigger than max_dt, 
        #we instead just jump to t+max_dt
        #and update the population at large.  If the step is smaller than max_dt
        #
        #()
        
        dt = random.expovariate(total_rate)
        if dt > max_dt:
            dt = max_dt
            t += dt
            deriv, FOI_Public = dPublicdt(Pop_status, (lambda_1, lambda_2, lambda_A, gamma_E, gamma_I1, gamma_I2, gamma_IA, q))  
            Pop_status = Pop_status+ dt *  deriv  
            lamb = (lambda_1*I1+lambda_2*I2)
            #FoHH = c_oHHx * lambda_1 * sum( sCohort.I1H for sCohort in subCohorts)
            #assert(FoHH==0)
            #deriv, FOI_oH = dCovC_HCW(oH_status, (c_H, lamb, K, c_HH, lambda_1, lambda_A, FoHH , gamma_E, gamma_I1, gamma_I2, gamma_IA, q))
            #oH_status =  oH_status + dt * deriv
            #for sCohort, PFOI, HFOI in zip(subCohorts, P_FOIs, HCW_FOIs):
            P_FOIs.append(myCohort.to_EP_rate/myCohort.SP)
            HCW_FOIs.append(myCohort.to_EH_rate/myCohort.SH)
            #CH, lamb, K, CHH, lambda1, FoHH , gammaE, gammaI1, gammaI2, gammaIA, q= parameters
            #Events: arrival, transmission, transition, departure
        else:
            t += dt
            deriv, FOI_Public = dPublicdt(Pop_status, (lambda_1, lambda_2, lambda_A, gamma_E, gamma_I1, gamma_I2, gamma_IA, q))  
            Pop_status = Pop_status + dt * deriv    
            #lamb = (lambda_1*I1+lambda_2*I2)
            # FoHH = c_oHHx * lambda_1 * sum( sCohort.I1H for sCohort in subCohorts)
            # deriv, FOI_oH = dCovC_HCW(oH_status, (c_H, lamb, K, c_HH, lambda_1, lambda_A, FoHH , gamma_E, gamma_I1, gamma_I2, gamma_IA, q))
            # oH_status +=  dt * deriv
            #CH, lamb, K, CHH, lambda1, FoHH , gammaE, gammaI1, gammaI2, gammaIA, q= parameters
            #for sCohort, PFOI, HFOI in zip(subCohorts, P_FOIs, HCW_FOIs):
            P_FOIs.append(myCohort.to_EP_rate/myCohort.SP)
            HCW_FOIs.append(myCohort.to_EH_rate/myCohort.SH)
            #r = random.random()*total_rate
            #for sCohort in subCohorts:
            #    r-= sCohort.total_rate
            #    if r<0:
            #        break
            myCohort.do_action()
        pop_FOI.append(FOI_Public)
#        oH_FOI.append(FOI_oH)
    #print('len', len(P_FOIs), len(HCW_FOIs))
    #print(len(subCohorts))
    return ts, pop_state, pop_FOI, Health_Facility_States, P_FOIs, HCW_FOIs
    
    
## MAIN CODE ##
if __name__ == '__main__':
    ###PARAMETERS###
    max_dt = 0.01  #maximum time step

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
    NumberPatients = 1000 #initial number of nonCovC patients
    NumberCohorts = 1  #number of nCovC subcohorts
#    Nhat = NumberPatients/NumberCohorts  #typical subCohort size
#    b = Nhat/ave_nonCOVIDduration   #entry rate required to have the typical stay be ave_nonCOVIDduration and size be Nhat.
    
    #oH_count = 0*1000  #initial number of CovC health care workers
    
    
    #oH_status = [oH_count, 0, 0, 0, 0, 0]  #S, E, I1, I2, IA, Rfor Covid-19 cohort of HCW
    
    Hx_count = int(round(NumberPatients/1))  # number of HCW in nonCov cohort Assuming 4 patients/HCW
    
    subCohorts = [subCohort(Hx_count/NumberCohorts,int(NumberPatients/NumberCohorts)) 
                for n in range(NumberCohorts)]
    Hx_statuses = [[int(Hx_count/NumberCohorts), 0, 0, 0, 0, 0] for n in range(NumberCohorts)]
    HCW_nCovC_trans_rates = [[] for n in range(NumberCohorts)]
    
    #P_statuses = Patient_nCovC_status = [np.array([int(NumberPatients/NumberCohorts), 0, 0]) for n in range(NumberCohorts)]
    Patient_nCovC_trans_rates = [[] for n in range(NumberCohorts)]
    
    Pop_status = np.array([0.99, 0.01, 0, 0, 0, 0])  #initial condition for population, normalized to 1]
    
    
    ##TRANSMISSION PARAMETERS##
    
    K =0*0.02  #transmission rate to COVID-19 HCWs from all infected patients combined (independent of number infected)
    
    c_P = 0.1# relative force of infection from general population I1 and IA to nCovC patients.
    c_H = 0.1# relative force of infection from general population infected to HCW.
    c_PP = 0.5# 
    c_HxP = 2
    c_PHx = 2 #relative force of infection from patients to HCWs in cohort compared to in general public
    c_HxHx = 1 # relative force of infection between HCWs in same subcohort
    c_HH = 0*0.05 #relative force of infection between HCWs within same cohort
    #c_oHHx = 0*0.01 #relative force of infection from CovC HCWs to nCovC HCWs (assume it's symmetric)
    assert(c_HH==0)
    omega = 0*0.05
    gamma_Q = 1./14
    ave_nonCOVIDduration = 14  #average hospital stay for non-COVID patient
    #print('fix line above')
    
    T=100 #end time
    
    ts, pop_state, pop_FOI, Health_Facility_States, P_FOIs, HCW_FOIs = simulate(
            lambda_1, lambda_2, lambda_A, omega, gamma_E, gamma_I1, gamma_I2, 
            gamma_IA, gamma_Q, q, 
            c_H, c_P,c_PP, c_HxP, c_PHx, c_HxHx, NumberPatients, #NumberCohorts, 
            T, ave_nonCOVIDduration=ave_nonCOVIDduration, PatientsPerHCW=4)
                
    #Patient_transmission_rates = [sum(L) for L in Patient_nCovC_trans_rates]
    #HCW_transmission_rates = [sum(L) for L in HCW_nCovC_trans_rates]
    plt.figure(3)
    plt.clf()
    nonCOVstates = [np.array(nonCOVstate) for nonCOVstate in Health_Facility_States]
    pop_state = np.array(pop_state)
    plt.plot(ts, pop_state[:,1], label = 'population E')
    plt.plot(ts, pop_state[:,0], label = 'population S')
    nonCOVstate =sum(nonCOVstates)
    
    #[sCohort.SP, sCohort.EP, sCohort.I1P, sCohort.IAP, sCohort.RP, sCohort.SH, sCohort.EH, sCohort.I1H, sCohort.IAH, sCohort.QH, sCohort.RH])
    
    plt.plot(ts, nonCOVstate[:,5]/nonCOVstate[0,5], '-.', label = 'S  HCW' )
    # plt.plot(ts, nonCOVstate[:,6]/nonCOVstate[0,5], '-.', label = 'E HCW')
    # plt.plot(ts, nonCOVstate[:,7]/nonCOVstate[0,5], '-.', label = 'I1 HCW')
    # plt.plot(ts, nonCOVstate[:,8]/nonCOVstate[0,5], '-.', label = 'IA HCW')
    # plt.plot(ts, nonCOVstate[:,9]/nonCOVstate[0,5], '-.', label = 'Q HCW')
    # plt.plot(ts, nonCOVstate[:,10]/nonCOVstate[0,5], '-.', label = 'R HCW')
    plt.plot(ts, nonCOVstate[:,0]/nonCOVstate[0,0], '--', label = 'S  Patient')
    plt.plot(ts, nonCOVstate[:,1]/nonCOVstate[0,0], '--', label = 'E  Patient')
    plt.plot(ts, nonCOVstate[:,2]/nonCOVstate[0,0], '--', label = 'I1  Patient')
    plt.plot(ts, nonCOVstate[:,3]/nonCOVstate[0,0], '--', label = 'IA  Patient')
    plt.plot(ts, nonCOVstate[:,4]/nonCOVstate[0,0], '--', label = 'R  Patient')
    plt.plot(ts, (nonCOVstate[:,0]+nonCOVstate[:,1]+nonCOVstate[:,2]+nonCOVstate[:,3]+nonCOVstate[:,4])/nonCOVstate[0,0])
    plt.legend()
    plt.show()  
    
    plt.figure(4)
    plt.clf()
    plt.plot(ts, pop_FOI, label = 'public FOI')
    for PFOI in P_FOIs:
        plt.plot(ts, PFOI, '--', label = 'subCohort Patient FOI')
    for HFOI in HCW_FOIs:
        plt.plot(ts, HFOI, '-.', label = 'subCohort HCW FOI')
    plt.legend()
    plt.show()
    #print(HCW_nCovC_status)