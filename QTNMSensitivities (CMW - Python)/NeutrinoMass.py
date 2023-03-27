#%% 

import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.stats import norm,halfnorm
import math
from math import log10,floor
import numba

#%% 

class DEopt: 
    
    def __init__(self,DE,vol,time,n,eff,b):
        
        self.DE = DE           # Analysis Window
        self.vol = vol         # Volume
        self.time = time       # Runtime (in years)
        self.n = n             # Number Density
        self.eff = eff         # Efficiency Factor
        self.b = b             # Background Rate
        
    def values(self):
        
        eta = 2.00e-13 #Branching Ratio in last eV
        self.b = 3.1536e+7*self.b # s^-1 to yr^-1 conversion
        self.n = self.n*(1.00e+6)
        tau_m = 12.3/np.log(2)
        R = (self.n*self.vol)/tau_m # Decay Rate
        r = self.eff*R*eta # Detected Decay Rate
        
        SD = []
        
        for i in range(len(self.DE)):
            term1 = r*self.time*self.DE[i]
            term2 = (self.b*self.time)/self.DE[i]
            factor = 2/(3*r*self.time)
            SD.append(factor*np.sqrt(term1+term2))
        
        return SD
          
    def min(self):
        
        eta = 2.00e-13 #Branching Ratio in last eV
        self.b = 3.1536e+7*self.b # s^-1 to yr^-1 conversion
        self.n = self.n*(1.00e+6)
        tau_m = 12.3/np.log(2)
        R = (self.n*self.vol)/tau_m # Decay Rate
        r = self.eff*R*eta # Detected Decay Rate
        
        SD = []
        
        for i in range(len(self.DE)):
            term1 = r*self.time*self.DE[i]
            term2 = (self.b*self.time)/self.DE[i]
            factor = 2/(3*r*self.time)
            SD.append(factor*np.sqrt(term1+term2))
        
        
        def rs(x,sig = 2):
            return round(x, sig-int(floor(log10(abs(x))))-1)
        
    
        print("DE_opt = (b/r)^(1/2) = ",rs(np.sqrt(self.b/r), sig = 5))
        
        for i in range(len(SD)):
            if SD[i] < SD[i+1]:
                print("Numerical DE_opt = ",rs((self.DE[i]), sig = 5))
                print("Standard Deviation for numerical DE_opt = ",rs(SD[i],sig = 5))
                break
            else:
                pass   
                
    def plot(self):
        
        DEvals = self.rangevals()
        fig = plt.figure(figsize = (10,8))
        plt.plot(self.DE,DEvals, 'b')
        plt.xlabel("Analysis Window $\\Delta E$", size = 18)
        plt.ylabel("Statistical Uncertainty $\\sigma_{m_\\beta^2}$(stat), eV$^{2}$", size = 18)
        plt.xscale("log")
        plt.yscale("log")
        plt.xticks(size = 16)
        plt.yticks(size = 16)
        plt.title(f"$\\sigma_{{m_\\beta{{^2}}}}$ against $\\Delta E$", size = 20)

class NMevars:
    
    def __init__(self, H3, n , vol , b , time , eff, B, rmsB , u_B , u_s, rmsI, u_I):
        
        # Key Experimental Parameters 
        
        self.H3,self.n,self.vol,self.B,self.time,self.eff,self.b = H3,n,vol,B,time,eff,b
        
        # Resolution in Systematic Uncertainties
        
        self.rmsB,self.rmsI = rmsB,rmsI
        
        # Uncertainties in Key Systematics.
        
        self.u_B,self.u_s,self.u_I = u_B,u_s,u_I
    
    def SDexposure(self):
    
        ''' Constants Relating to Relativistic Dynamics of Decay Electrons'''
        
        # Mass, Cyclotron Frequency, Speed of Light, Tritium Endpoint Energy
        
        mass,f_c,c,E_0 = 9.109e-31,26.5e+9,2.99792458e8,18.6*(1.60e-16) 
        
        # Electron Rest Mass Energy, Lorentz Factor, Relativistic Ratio          
        
        Emc= mass*(c**2)                        
        gamma = 1 + E_0/Emc
        beta = np.sqrt(1-1/(gamma**2))
        ''' -------------------- Tritium Decay Dynamics ---------------------'''
        
        self.n = self.n*(1.00e+6)             # Number Density of Tritium
        yr,eta = 3.1536e+7,2.00e-13
        tau_m = 12.3/np.log(2) # Branching Ratio in Last eV, Mean Lifetime of Tritium
        
        rate = []
        decay_factor = self.eff*(self.n/tau_m)*eta
        
        self.b = yr*self.b
        
        for i in range(len(self.vol)):
            rate.append(decay_factor*self.vol[i]) # Detected Decay Rate due to Efficiency (or "Solid Pitch Angle")
        
        ''' -------------------- Statistical Uncertainty ----------------------'''
        
        ustat,factor,term1,term2 = [],[],[],[]
        DE = 1.0
        #for i in range(len(self.vol)):
            #DE.append(np.sqrt(self.b/rate[i]))
        
        for j in range(len(rate)):
            factor.append(2/(3*rate[j]*self.time))
            term1.append((rate[j]*self.time*DE))
            term2.append((self.b*self.time)/DE)
            ustat.append(factor[j]*np.sqrt(term1[j]+term2[j]))
        

        ''' -------------------- Systematic Uncertainties ----------------------'''
        
        E_end =  18.6e+3                # Endpoint Energy(keV),
        Ey = E_end/(gamma-1)            # Energy Factor
        instr_val = 1.0                 #Instrument Resolution
        instr_res = self.rmsI*instr_val
        ME,IE = Ey*(self.rmsB/self.B), Ey*(instr_res/instr_val) # Magnetic Field and Instrument Systematic Error        
        
        u1,u2,u3 = self.u_B , self.u_s , self.u_I               # Percentage Uncertainties (0.00 - 1.00)

        if self.H3 == "A":
            
            cs = 9.0*(1.00e-19)
            sigma0 = cs*(1.00e-4)
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)
        
            systematics = [ME,SE,IE]
            uncertainty = [u1,u2,u3]
            usyst = []
            for i in range(len(systematics)):
                usyst.append((uncertainty[i]**2)*(systematics[i]**4))
            
            usyst = 4*np.sqrt(np.sum((usyst)))
            
            total = []                          # Total Uncertainty
            for j in range(len(ustat)):
                total.append(np.sqrt(ustat[j]**2 + usyst**2))
            
            return total
        
        if self.H3 == "M": 
            
            cs = 3.4*(1.00e-18)
            sigma0 = cs*(1.00e-4)
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)
            FSS = 0.30 #0.2978723404255319
                    
            systematics = [ME,SE,FSS,IE]
            uncertainty = [u1,u2,0.01,u3]
            usyst = []
            
            for i in range(len(systematics)):
                usyst.append((uncertainty[i]**2)*(systematics[i]**4))
            
            usyst = 4*np.sqrt(np.sum((usyst)))
            
            total = []                          # Total Uncertainty
            for j in range(len(ustat)):
                total.append(np.sqrt(ustat[j]**2 + usyst**2))
            
            return total
            
    def MSexposure(self,confidence): 
        valsSD = self.SDexposure()
        self.confidence = confidence
        expectation = 0
        # Function will depend on standard deviation
        msv= [] # Mass squared values
        for i in range(len(valsSD)):
            cl = halfnorm.interval(confidence,expectation,valsSD[i])
            msv.append(cl[1])
        mass_vals = np.sqrt(msv)
        return mass_vals
        
    def SDndensity(self):
        ''' Constants Relating to Relativistic Dynamics of Decay Electrons'''
        
        # Mass, Cyclotron Frequency, Speed of Light, Tritium Endpoint Energy
        
        mass,f_c,c,E_0 = 9.109e-31,26.5e+9,2.99792458e8,18.6*(1.60e-16) 
        
        # Electron Rest Mass Energy, Lorentz Factor, Relativistic Ratio          
        
        Emc= mass*(c**2)                        
        gamma = 1 + E_0/Emc
        beta = np.sqrt(1-1/(gamma**2))
        ''' -------------------- Tritium Decay Dynamics ---------------------'''
        
        yr,eta = 3.1536e+7,2.00e-13
        tau_m = 12.3/np.log(2) # Branching Ratio in Last eV, Mean Lifetime of Tritium
        
        rate = []
        decay_factor = self.eff*(self.vol/tau_m)*eta
        
        self.b = yr*self.b
        nm3 = [] # Number Density in Cubic Meters

        for i in range(len(self.n)):
            nm3.append(self.n[i]*1.00e+6)
        for i in range(len(self.n)):
            rate.append(decay_factor*nm3[i]) # Detected Decay Rate due to Efficiency (or "Solid Pitch Angle")
        
        ''' -------------------- Statistical Uncertainty ----------------------'''
        
        ustat,factor,term1,term2 = [],[],[],[]
        DE = 1.0
        
        for j in range(len(self.n)):
            factor.append(2/(3*rate[j]*self.time))
            term1.append((rate[j]*self.time*DE))
            term2.append((self.b*self.time)/DE)
            ustat.append(factor[j]*np.sqrt(term1[j]+term2[j]))
            
        ''' -------------------- Systematic Uncertainties ----------------------'''
        
        E_end =  18.6e+3                # Endpoint Energy(keV),
        Ey = E_end/(gamma-1)            # Energy Factor
        instr_val = 1.0                 #Instrument Resolution
        instr_res = self.rmsI*instr_val
        ME,IE = Ey*(self.rmsB/self.B), Ey*(instr_res/instr_val) # Magnetic Field and Instrument Systematic Error        
        
        u1,u2,u3 = self.u_B , self.u_s , self.u_I               # Percentage Uncertainties (0.00 - 1.00)

        if self.H3 == "A":
            
            scat2 = 2*np.pi*f_c
            cs = 9.0*(1.00e-19)
            sigma0 = cs*(1.00e-4)
            
            SE = []
            usyst1,usyst2 = [], []
            scat1 = []
            for i in range(len(self.n)):
                scat1.append(beta*c*sigma0*nm3[i])
            for i in range(len(self.n)):
                SE.append(Ey*(scat1[i]/scat2))
                usyst1.append((ME**4)*(u1**2) + (SE[i]**4)*(u2**2) + (IE**4)*(u3**2))
                usyst2.append(4*np.sqrt(usyst1[i]))
            
            total = []                          # Total Uncertainty
            for j in range(len(self.n)):
                total.append(np.sqrt(ustat[j]**2 + usyst2[j]**2))
            
            return total
        
        if self.H3 == "M":
            
            
            scat2 = 2*np.pi*f_c
            cs = 3.4*(1.00e-18)
            sigma0 = cs*(1.00e-4)
            FSS = 0.30
            
            SE = []
            usyst1,usyst2 = [], []
            scat1 = []
            for i in range(len(self.n)):
                scat1.append(beta*c*sigma0*nm3[i])
            for i in range(len(self.n)):
                SE.append(Ey*(scat1[i]/scat2))
                usyst1.append((ME**4)*(u1**2) + (SE[i]**4)*(u2**2) + (IE**4)*(u3**2) + (FSS**4)*(0.01**2))
                usyst2.append(4*np.sqrt(usyst1[i]))
            
            total = []                          # Total Uncertainty
            for j in range(len(self.n)):
                total.append(np.sqrt(ustat[j]**2 + usyst2[j]**2))
            
            return total
               
    def MSndensity(self,confidence):
        valsSD = self.SDndensity()
        self.confidence = confidence
        expectation = 0
        # Function will depend on standard deviation
        msv= [] # Mass squared values
        for i in range(len(valsSD)):
            cl = halfnorm.interval(confidence,expectation,valsSD[i])
            msv.append(cl[1])
        mass_vals = np.sqrt(msv)
        return mass_vals

    def MSndensity_opt(self,confidence):
        self.confidence = confidence
        valsMS = self.MSndensity(self.confidence)
        expectation = 0
        # Function will depend on standard deviation
        msv= [] # Mass squared values
       
        for i in range(len(valsMS)):
            if valsMS[i] < valsMS[i+1]:
                print("Optimum Number Density", self.n[i])
                print("Minimum Neutrino Mass at Optimum n ",valsMS[i], "eV")
                break
         
    def SDbackground(self):
        
        ''' Constants Relating to Relativistic Dynamics of Decay Electrons'''
        
        # Mass, Cyclotron Frequency, Speed of Light, Tritium Endpoint Energy
        
        mass,f_c,c,E_0 = 9.109e-31,26.5e+9,2.99792458e8,18.6*(1.60e-16) 
        
        # Electron Rest Mass Energy, Lorentz Factor, Relativistic Ratio          
        
        Emc= mass*(c**2)                        
        gamma = 1 + E_0/Emc
        beta = np.sqrt(1-1/(gamma**2))
        ''' -------------------- Tritium Decay Dynamics ---------------------'''
        
        self.n = self.n*(1.00e+6)             # Number Density of Tritium
        yr,eta = 3.1536e+7,2.00e-13
        tau_m = 12.3/np.log(2) # Branching Ratio in Last eV, Mean Lifetime of Tritium
    
        rate = self.eff*(self.n*self.vol/tau_m)*eta
        
        back = []
        for i in range(len(self.b)):
            back.append(yr*self.b[i])
            
        
    
        ''' -------------------- Statistical Uncertainty ----------------------'''
        
        DE = 1.0
        factor = 2/(3*rate*self.time)
        term1 = rate*self.time*DE
        ustat,term2 = [],[]
        
        for j in range(len(back)):
            term2.append((back[j]*self.time)/DE)
            ustat.append(factor*np.sqrt(term1+term2[j]))
        
        ''' -------------------- Systematic Uncertainties ----------------------'''
        
        E_end =  18.6e+3                # Endpoint Energy(keV),
        Ey = E_end/(gamma-1)            # Energy Factor
        instr_val = 1.0                 #Instrument Resolution
        instr_res = self.rmsI*instr_val
        ME,IE = Ey*(self.rmsB/self.B), Ey*(instr_res/instr_val) # Magnetic Field and Instrument Systematic Error        
        
        u1,u2,u3 = self.u_B , self.u_s , self.u_I               # Percentage Uncertainties (0.00 - 1.00)

        if self.H3 == "A":
            
            cs = 9.0*(1.00e-19)
            sigma0 = cs*(1.00e-4)
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)
        
            systematics = [ME,SE,IE]
            uncertainty = [u1,u2,u3]
            usyst = []
            for i in range(len(systematics)):
                usyst.append((uncertainty[i]**2)*(systematics[i]**4))
            
            usyst = 4*np.sqrt(np.sum((usyst)))
            
            total = []                          # Total Uncertainty
            for j in range(len(ustat)):
                total.append(np.sqrt(ustat[j]**2 + usyst**2))
            
            return total
        
        if self.H3 == "M": 
            
            cs = 3.4*(1.00e-18)
            sigma0 = cs*(1.00e-4)
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)
            FSS = 0.30 #0.2978723404255319
                  
            systematics = [ME,SE,FSS,IE]
            uncertainty = [u1,u2,0.01,u3]
            usyst = []
            
            for i in range(len(systematics)):
                usyst.append((uncertainty[i]**2)*(systematics[i]**4))
            
            usyst = 4*np.sqrt(np.sum((usyst)))
            
            total = []                          # Total Uncertainty
            for j in range(len(ustat)):
                total.append(np.sqrt(ustat[j]**2 + usyst**2))
            
            return total
        
    def MSbackground(self,confidence):
        valsSD = self.SDbackground()
        self.confidence = confidence
        expectation = 0
        # Function will depend on standard deviation
        msv= [] # Mass squared values
        for i in range(len(valsSD)):
            cl = halfnorm.interval(confidence,expectation,valsSD[i])
            msv.append(cl[1])
        mass_vals = np.sqrt(msv)
        mass_vals1 = []
        for i in range(len(mass_vals)):
            mass_vals1.append(mass_vals[i]*1000)
        return mass_vals1
    
class NMuvars1:
    
    def __init__(self, H3, n , vol , b , time , eff, B, rmsB , u_B , u_s, rmsI, u_I):
        
        # Experimental Parameters
        self.H3, self.n , self.vol, self.time, self.eff , self.B,self.b = H3 , n , vol, time, eff,B,b
        
        # Variation (Uncertainty) with the instrumentals

        self.rmsB = rmsB
        self.rmsI = rmsI

        # Uncertainty on rms Variation 

        self.u_B = u_B
        self.u_I = u_I
        self.u_s = u_s
        
    def SDscattering(self):
       
        ''' Constants Relating to Relativistic Dynamics of Decay Electrons'''
        
        # Mass, Cyclotron Frequency, Speed of Light, Tritium Endpoint Energy
        
        mass,f_c,c,E_0 = 9.109e-31,26.5e+9,2.99792458e8,18.6*(1.60e-16) 
        
        # Electron Rest Mass Energy, Lorentz Factor, Relativistic Ratio          
        
        Emc= mass*(c**2)                        
        gamma = 1 + E_0/Emc
        beta = np.sqrt(1-1/(gamma**2))
        ''' -------------------- Tritium Decay Dynamics ---------------------'''
        
        yr,eta = 3.1536e+7,2.00e-13
        tau_m = 12.3/np.log(2) # Branching Ratio in Last eV, Mean Lifetime of Tritium
        
        self.b = yr*self.b
        self.n = self.n*(1.00e+6)

        rate = self.eff*(self.vol*self.n/tau_m)*eta
        
        ''' -------------------- Statistical Uncertainty ----------------------'''
        
        DE = 1.00e0
        factor = 2/(3*rate*self.time)
        term1 = rate*self.time*DE
        term2 = (self.b*self.time)/DE
        ustat = factor*np.sqrt(term1 + term2)
       
        ''' -------------------- Systematic Uncertainties ----------------------'''
        
        E_end =  18.6e+3                # Endpoint Energy(keV),
        Ey = E_end/(gamma-1)            # Energy Factor
        instr_val = 1.0                 #Instrument Resolution
        instr_res = self.rmsI*instr_val
        ME,IE = Ey*(self.rmsB/self.B), Ey*(instr_res/instr_val) # Magnetic Field and Instrument Systematic Error        
        
        if self.H3 == "A":
            
            cs = 9.0*(1.00e-19)
            sigma0 = cs*(1.00e-4)
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)

            SEsyst = []
            for i in range(len(self.u_s)):
                SEsyst.append((SE**4)*(self.u_s[i]**2))

            usyst1,usyst2 = [],[]
         
            for j in range(len(self.u_s)):
                usyst1.append((ME**4)*(self.u_B**2) + SEsyst[j] + (IE**4)*(self.u_I**2))
                usyst2.append(4*np.sqrt(usyst1[j]))
            
            total = []                          # Total Uncertainty
            for j in range(len(self.u_s)):
                total.append(np.sqrt(ustat**2 + usyst2[j]**2))
            
            return total
        
        if self.H3 == "M":
            
            scat2 = 2*np.pi*f_c
            cs = 3.4*(1.00e-18)
            sigma0 = cs*(1.00e-4)
            FSS = 0.30
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)

            SEsyst = []
            for i in range(len(self.u_s)):
                SEsyst.append((SE**4)*(self.u_s[i]**2))

            usyst1,usyst2 = [],[]
         
            for j in range(len(self.u_s)):
                usyst1.append((ME**4)*(self.u_B**2) + SEsyst[j] + (IE**4)*(self.u_I**2)+(FSS**4)*(0.01**2))
                usyst2.append(4*np.sqrt(usyst1[j]))
            
            total = []                          # Total Uncertainty
            for j in range(len(self.u_s)):
                total.append(np.sqrt(ustat**2 + usyst2[j]**2))
            
            return total
    
    def MSscattering(self,confidence):
        valsSD = self.SDscattering()
        self.confidence = confidence
        expectation = 0
        # Function will depend on standard deviation
        msv= [] # Mass squared values
        for i in range(len(valsSD)):
            cl = halfnorm.interval(confidence,expectation,valsSD[i])
            msv.append(cl[1])
        mass_vals = np.sqrt(msv)
        return mass_vals
    
    def SDmagnetic(self):
        
        ''' Constants Relating to Relativistic Dynamics of Decay Electrons'''
        
        # Mass, Cyclotron Frequency, Speed of Light, Tritium Endpoint Energy
        
        mass,f_c,c,E_0 = 9.109e-31,26.5e+9,2.99792458e8,18.6*(1.60e-16) 
        
        # Electron Rest Mass Energy, Lorentz Factor, Relativistic Ratio          
        
        Emc= mass*(c**2)                        
        gamma = 1 + E_0/Emc
        beta = np.sqrt(1-1/(gamma**2))
        ''' -------------------- Tritium Decay Dynamics ---------------------'''
        
        yr,eta = 3.1536e+7,2.00e-13
        tau_m = 12.3/np.log(2) # Branching Ratio in Last eV, Mean Lifetime of Tritium
        
        self.b = yr*self.b
        self.n = self.n*(1.00e+6)

        rate = self.eff*(self.vol*self.n/tau_m)*eta
        
        ''' -------------------- Statistical Uncertainty ----------------------'''
        
        DE = 1.00e0
        factor = 2/(3*rate*self.time)
        term1 = rate*self.time*DE
        term2 = (self.b*self.time)/DE
        ustat = factor*np.sqrt(term1 + term2)
       
        ''' -------------------- Systematic Uncertainties ----------------------'''
        
        E_end =  18.6e+3                # Endpoint Energy(keV),
        Ey = E_end/(gamma-1)            # Energy Factor
        instr_val = 1.0                 #Instrument Resolution
        instr_res = self.rmsI*instr_val
        ME,IE = Ey*(self.rmsB/self.B), Ey*(instr_res/instr_val) # Magnetic Field and Instrument Systematic Error  
        
        if self.H3 == "A":
            
            cs = 9.0*(1.00e-19)
            sigma0 = cs*(1.00e-4)
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)

            MEsyst = []
            for i in range(len(self.u_B)):
                MEsyst.append((ME**4)*(self.u_B[i]**2))

            usyst1,usyst2 = [],[]
         
            for j in range(len(self.u_B)):
                usyst1.append((SE**4)*(self.u_s**2) + MEsyst[j] + (IE**4)*(self.u_I**2))
                usyst2.append(4*np.sqrt(usyst1[j]))
            
            total = []                          # Total Uncertainty
            for j in range(len(self.u_B)):
                total.append(np.sqrt(ustat**2 + usyst2[j]**2))
            
            return total
        
        if self.H3 == "M":
            
            scat2 = 2*np.pi*f_c
            cs = 3.4*(1.00e-18)
            sigma0 = cs*(1.00e-4)
            FSS = 0.30
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)

            MEsyst = []
            for i in range(len(self.u_B)):
                MEsyst.append((ME**4)*(self.u_B[i]**2))

            usyst1,usyst2 = [],[]
         
            for j in range(len(self.u_B)):
                usyst1.append((SE**4)*(self.u_s**2) + MEsyst[j] + (IE**4)*(self.u_I**2)+(FSS**4)*(0.01**2))
                usyst2.append(4*np.sqrt(usyst1[j]))
            
            total = []                          # Total Uncertainty
            for j in range(len(self.u_B)):
                total.append(np.sqrt(ustat**2 + usyst2[j]**2))
            
            return total
    
    def MSmagnetic(self,confidence):
        valsSD = self.SDmagnetic()
        self.confidence = confidence
        expectation = 0
        # Function will depend on standard deviation
        msv= [] # Mass squared values
        for i in range(len(valsSD)):
            cl = halfnorm.interval(confidence,expectation,valsSD[i])
            msv.append(cl[1])
        mass_vals = np.sqrt(msv)
        return mass_vals
    
    def SDinstrument(self):
        ''' Constants Relating to Relativistic Dynamics of Decay Electrons'''
        
        # Mass, Cyclotron Frequency, Speed of Light, Tritium Endpoint Energy
        
        mass,f_c,c,E_0 = 9.109e-31,26.5e+9,2.99792458e8,18.6*(1.60e-16) 
        
        # Electron Rest Mass Energy, Lorentz Factor, Relativistic Ratio          
        
        Emc= mass*(c**2)                        
        gamma = 1 + E_0/Emc
        beta = np.sqrt(1-1/(gamma**2))
        ''' -------------------- Tritium Decay Dynamics ---------------------'''
        
        yr,eta = 3.1536e+7,2.00e-13
        tau_m = 12.3/np.log(2) # Branching Ratio in Last eV, Mean Lifetime of Tritium
        
        self.b = yr*self.b
        self.n = self.n*(1.00e+6)

        rate = self.eff*(self.vol*self.n/tau_m)*eta
        
        ''' -------------------- Statistical Uncertainty ----------------------'''
        
        DE = 1.00e0
        factor = 2/(3*rate*self.time)
        term1 = rate*self.time*DE
        term2 = (self.b*self.time)/DE
        ustat = factor*np.sqrt(term1 + term2)
       
        ''' -------------------- Systematic Uncertainties ----------------------'''
        
        E_end =  18.6e+3                # Endpoint Energy(keV),
        Ey = E_end/(gamma-1)            # Energy Factor
        instr_val = 1.0                 #Instrument Resolution
        instr_res = self.rmsI*instr_val
        ME,IE = Ey*(self.rmsB/self.B), Ey*(instr_res/instr_val) # Magnetic Field and Instrument Systematic Error  
        
        if self.H3 == "A":
            
            cs = 9.0*(1.00e-19)
            sigma0 = cs*(1.00e-4)
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)

            IEsyst = []
            for i in range(len(self.u_I)):
                IEsyst.append((IE**4)*(self.u_I[i]**2))

            usyst1,usyst2 = [],[]
         
            for j in range(len(self.u_I)):
                usyst1.append((SE**4)*(self.u_s**2) + IEsyst[j] + (ME**4)*(self.u_B**2))
                usyst2.append(4*np.sqrt(usyst1[j]))
            
            total = []                          # Total Uncertainty
            for j in range(len(self.u_I)):
                total.append(np.sqrt(ustat**2 + usyst2[j]**2))
            
            return total
        
        if self.H3 == "M":
            
            scat2 = 2*np.pi*f_c
            cs = 3.4*(1.00e-18)
            sigma0 = cs*(1.00e-4)
            FSS = 0.30
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)

            IEsyst = []
            for i in range(len(self.u_I)):
                IEsyst.append((IE**4)*(self.u_I[i]**2))

            usyst1,usyst2 = [],[]
         
            for j in range(len(self.u_I)):
                usyst1.append((SE**4)*(self.u_s**2) + IEsyst[j] + (ME**4)*(self.u_B**2)+(FSS**4)*(0.01**2))
                usyst2.append(4*np.sqrt(usyst1[j]))
            
            total = []                          # Total Uncertainty
            for j in range(len(self.u_I)):
                total.append(np.sqrt(ustat**2 + usyst2[j]**2))
            
            return total
    
    def MSinstrument(self,confidence):
        valsSD = self.SDinstrument()
        self.confidence = confidence
        expectation = 0
        # Function will depend on standard deviation
        msv= [] # Mass squared values
        for i in range(len(valsSD)):
            cl = halfnorm.interval(confidence,expectation,valsSD[i])
            msv.append(cl[1])
        mass_vals = np.sqrt(msv)
        return mass_vals
 
class NMuvars2:
    
    def __init__(self, H3, n , vol , b , time , eff, B, rmsB , u_B , u_s, rmsI, u_I):
        
        # Experimental Parameters
        self.H3, self.n , self.vol, self.time, self.eff , self.B,self.b = H3 , n , vol, time, eff,B,b
        
        # Variation (Uncertainty) with the instrumentals

        self.rmsB = rmsB
        self.rmsI = rmsI

        # Uncertainty on rms Variation 

        self.u_B = u_B
        self.u_I = u_I
        self.u_s = u_s
        
    def SDres_magnetic(self):
        
        ''' Constants Relating to Relativistic Dynamics of Decay Electrons'''
        
        # Mass, Cyclotron Frequency, Speed of Light, Tritium Endpoint Energy
        
        mass,f_c,c,E_0 = 9.109e-31,26.5e+9,2.99792458e8,18.6*(1.60e-16) 
        
        # Electron Rest Mass Energy, Lorentz Factor, Relativistic Ratio          
        
        Emc= mass*(c**2)                        
        gamma = 1 + E_0/Emc
        beta = np.sqrt(1-1/(gamma**2))
        ''' -------------------- Tritium Decay Dynamics ---------------------'''
        
        yr,eta = 3.1536e+7,2.00e-13
        tau_m = 12.3/np.log(2) # Branching Ratio in Last eV, Mean Lifetime of Tritium
        
        self.b = yr*self.b
        self.n = self.n*(1.00e+6)

        rate = self.eff*(self.vol*self.n/tau_m)*eta
        
        ''' -------------------- Statistical Uncertainty ----------------------'''
        
        DE = 1.00e0
        factor = 2/(3*rate*self.time)
        term1 = rate*self.time*DE
        term2 = (self.b*self.time)/DE
        ustat = factor*np.sqrt(term1 + term2)
       
        ''' -------------------- Systematic Uncertainties ----------------------'''
        
        E_end =  18.6e+3                # Endpoint Energy(keV),
        Ey = E_end/(gamma-1)            # Energy Factor
        instr_val = 1.0                 #Instrument Resolution
        instr_res = self.rmsI*instr_val
        
        ME = []
        for i in range(len(self.rmsB)):
            ME.append(Ey*(self.rmsB[i]/self.B))
            
        IE = Ey*(instr_res/instr_val) # Magnetic Field and Instrument Systematic Error  
        
        if self.H3 == "A":
            
            cs = 9.0*(1.00e-19)
            sigma0 = cs*(1.00e-4)
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)
            MEsyst = []
            for i in range(len(self.rmsB)):
                MEsyst.append((ME[i]**4)*(self.u_B**2))

            usyst1,usyst2 = [],[]
         
            for j in range(len(self.rmsB)):
                usyst1.append((SE**4)*(self.u_s**2) + MEsyst[j] + (IE**4)*(self.u_I**2))
                usyst2.append(4*np.sqrt(usyst1[j]))
            
            total = []                          # Total Uncertainty
            for j in range(len(self.rmsB)):
                total.append(np.sqrt(ustat**2 + usyst2[j]**2))
            
            return total
        
        if self.H3 == "M":
            
            scat2 = 2*np.pi*f_c
            cs = 3.4*(1.00e-18)
            sigma0 = cs*(1.00e-4)
            FSS = 0.30
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)

            MEsyst = []
            for i in range(len(self.rmsB)):
                MEsyst.append((ME[i]**4)*(self.u_B**2))

            usyst1,usyst2 = [],[]
         
            for j in range(len(self.rmsB)):
                usyst1.append((SE**4)*(self.u_s**2) + MEsyst[j] + (IE**4)*(self.u_I**2)+(FSS**4)*(0.01**2))
                usyst2.append(4*np.sqrt(usyst1[j]))
            
            total = []                          # Total Uncertainty
            for j in range(len(self.u_B)):
                total.append(np.sqrt(ustat**2 + usyst2[j]**2))
            
            return total
    
    def MSres_magnetic(self,confidence):
        valsSD = self.SDres_magnetic()
        self.confidence = confidence
        expectation = 0
        # Function will depend on standard deviation
        msv= [] # Mass squared values
        for i in range(len(valsSD)):
            cl = halfnorm.interval(confidence,expectation,valsSD[i])
            msv.append(cl[1])
        mass_vals = np.sqrt(msv)
        return mass_vals
    
    def SDres_instrument(self):
        ''' Constants Relating to Relativistic Dynamics of Decay Electrons'''
        
        # Mass, Cyclotron Frequency, Speed of Light, Tritium Endpoint Energy
        
        mass,f_c,c,E_0 = 9.109e-31,26.5e+9,2.99792458e8,18.6*(1.60e-16) 
        
        # Electron Rest Mass Energy, Lorentz Factor, Relativistic Ratio          
        
        Emc= mass*(c**2)                        
        gamma = 1 + E_0/Emc
        beta = np.sqrt(1-1/(gamma**2))
        ''' -------------------- Tritium Decay Dynamics ---------------------'''
        
        yr,eta = 3.1536e+7,2.00e-13
        tau_m = 12.3/np.log(2) # Branching Ratio in Last eV, Mean Lifetime of Tritium
        
        self.b = yr*self.b
        self.n = self.n*(1.00e+6)

        rate = self.eff*(self.vol*self.n/tau_m)*eta
        
        ''' -------------------- Statistical Uncertainty ----------------------'''
        
        DE = 1.00e0
        factor = 2/(3*rate*self.time)
        term1 = rate*self.time*DE
        term2 = (self.b*self.time)/DE
        ustat = factor*np.sqrt(term1 + term2)
       
        ''' -------------------- Systematic Uncertainties ----------------------'''
        
        E_end =  18.6e+3                # Endpoint Energy(keV),
        Ey = E_end/(gamma-1)            # Energy Factor
        instr_val = 1.0                 #Instrument Resolution
        #instr_res = self.rmsI*instr_val
        ME = Ey*(self.rmsB/self.B)
        IE = []
        for i in range(len(self.rmsI)):
            IE.append(Ey*(self.rmsI[i]*instr_val))
        
        if self.H3 == "A":
            
            cs = 9.0*(1.00e-19)
            sigma0 = cs*(1.00e-4)
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)
            IEsyst = []
            for i in range(len(self.rmsI)):
                IEsyst.append((IE[i]**4)*(self.u_I**2))

            usyst1,usyst2 = [],[]
         
            for j in range(len(self.rmsI)):
                usyst1.append((SE**4)*(self.u_s**2) + IEsyst[j] + (ME**4)*(self.u_B**2))
                usyst2.append(4*np.sqrt(usyst1[j]))
            
            total = []                          # Total Uncertainty
            for j in range(len(self.rmsI)):
                total.append(np.sqrt(ustat**2 + usyst2[j]**2))
            
            return total
        
        if self.H3 == "M":
            
            FSS = 0.30
            cs = 9.0*(1.00e-19)
            sigma0 = cs*(1.00e-4)
            scat1 = beta*c*sigma0*self.n
            scat2 = 2*np.pi*f_c
            SE = Ey*(scat1/scat2)
            
            IEsyst = []
            for i in range(len(self.rmsI)):
                IEsyst.append((IE[i]**4)*(self.u_I**2))
           
            usyst1,usyst2 = [],[]
         
            for j in range(len(self.rmsI)):
                usyst1.append((SE**4)*(self.u_s**2) + IEsyst[j] + (ME**4)*(self.u_B**2)++(FSS**4)*(0.01**2))
                usyst2.append(4*np.sqrt(usyst1[j]))
            
            total = []                          # Total Uncertainty
            for j in range(len(self.rmsI)):
                total.append(np.sqrt(ustat**2 + usyst2[j]**2))
             
            return total
    
    def MSres_instrument(self,confidence):
        valsSD = self.SDres_instrument()
        self.confidence = confidence
        expectation = 0
        # Function will depend on standard deviation
        msv= [] # Mass squared values
        for i in range(len(valsSD)):
            cl = halfnorm.interval(confidence,expectation,valsSD[i])
            msv.append(cl[1])
        mass_vals = np.sqrt(msv)
        return mass_vals

def NMfinmass(H3, n , vol , b , time , eff, B, rmsB , u_B , u_s, rmsI, u_I,confidence):
    
    ''' Constants Relating to Relativistic Dynamics of Decay Electrons'''
        
    # Mass, Cyclotron Frequency, Speed of Light, Tritium Endpoint Energy
    
    mass,f_c,c,E_0 = 9.109e-31,26.5e+9,2.99792458e8,18.6*(1.60e-16) 
    
    # Electron Rest Mass Energy, Lorentz Factor, Relativistic Ratio          
        
    Emc= mass*(c**2)                        
    gamma = 1 + E_0/Emc
    beta = np.sqrt(1-1/(gamma**2))
    ''' -------------------- Tritium Decay Dynamics ---------------------'''
    
    yr,eta = 3.1536e+7,2.00e-13
    tau_m = 12.3/np.log(2) # Branching Ratio in Last eV, Mean Lifetime of Tritium
    
    b = yr*b
    n = n*(1.00e+6)

    rate = eff*(vol*n/tau_m)*eta
    
    ''' -------------------- Statistical Uncertainty ----------------------'''
    
    DE = 1.00e0
    factor = 2/(3*rate*time)
    term1 = rate*time*DE
    term2 = (b*time)/DE
    ustat = factor*np.sqrt(term1 + term2)
    
    ''' -------------------- Systematic Uncertainties ----------------------'''
    
    E_end =  18.6e+3                # Endpoint Energy(keV),
    Ey = E_end/(gamma-1)            # Energy Factor
    instr_val = 1.0                 #Instrument Resolution
    instr_res = rmsI*instr_val
    ME,IE = Ey*(rmsB/B), Ey*(instr_res/instr_val) # Magnetic Field and Instrument Systematic Error        
    
    u1,u2,u3 = u_B , u_s , u_I               # Percentage Uncertainties (0.00 - 1.00)

    if H3 == "A":
        
        cs = 9.0*(1.00e-19)
        sigma0 = cs*(1.00e-4)
        scat1 = beta*c*sigma0*n
        scat2 = 2*np.pi*f_c
        SE = Ey*(scat1/scat2)
    
        systematics = [ME,SE,IE]
        uncertainty = [u1,u2,u3]
        usyst = []
        for i in range(len(systematics)):
            usyst.append((uncertainty[i]**2)*(systematics[i]**4))
        
        usyst = 4*np.sqrt(np.sum((usyst)))
        total = np.sqrt(ustat**2 + usyst**2)
        expectation = 0
        cl = halfnorm.interval(confidence,expectation,total)
        mass_vals = np.sqrt(cl[1])
        print(mass_vals*1000,"meV")
    if H3 == "M": 
        
        cs = 3.4*(1.00e-18)
        sigma0 = cs*(1.00e-4)
        scat1 = beta*c*sigma0*n
        scat2 = 2*np.pi*f_c
        SE = Ey*(scat1/scat2)
        FSS = 0.30 #0.2978723404255319
                
        systematics = [ME,SE,FSS,IE]
        uncertainty = [u1,u2,0.01,u3]
        usyst = []
        
        for i in range(len(systematics)):
            usyst.append((uncertainty[i]**2)*(systematics[i]**4))
        
        usyst = 4*np.sqrt(np.sum((usyst)))
        total = np.sqrt(ustat**2 + usyst**2)
        expectation = 0
        cl = halfnorm.interval(confidence,expectation,total)
        mass_vals = np.sqrt(cl[1])
        print(mass_vals*1000,"meV")
    
    
# %%
