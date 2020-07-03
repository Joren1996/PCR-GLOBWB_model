# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 16:07:09 2020

@author: joren

Seperate Code file with the Richard's model, unsaturated zone temperature and the energy balance included.

"""

import numpy as np
import pcraster as pcr
from scipy.integrate import odeint
from scipy import optimize
import multiprocessing as mp

import logging
logger = logging.getLogger(__name__)


def initialiseRichards(self):
    msg='Begin initialisation Richards model!'
    logger.info(msg)
    
    #Convert relevant variables to numpy.
    global storCapUpp
    kSatUpp_numpy=pcr.pcr2numpy(self.parameters.kSatUpp, np.nan)
    kSatLow_numpy=pcr.pcr2numpy(self.parameters.kSatLow, np.nan)
    satVolMoistContUpp_numpy=pcr.pcr2numpy(self.parameters.satVolMoistContUpp, np.nan)
    satVolMoistContLow_numpy=pcr.pcr2numpy(self.parameters.satVolMoistContLow, np.nan)
    resVolMoistContUpp_numpy=pcr.pcr2numpy(self.parameters.resVolMoistContUpp, np.nan)
    resVolMoistContLow_numpy=pcr.pcr2numpy(self.parameters.resVolMoistContLow, np.nan)
    thickUpp=pcr.pcr2numpy(self.parameters.thickUpp, np.nan)
    thickLow=pcr.pcr2numpy(self.parameters.thickLow, np.nan)
    storCapUpp=pcr.pcr2numpy(self.parameters.storCapUpp, np.nan)
    storUpp_numpy=pcr.pcr2numpy(self.storUpp, np.nan)
    storLow_numpy=pcr.pcr2numpy(self.storLow, np.nan)

    #Dimensions of the original basin
    self.basin_shape=np.shape(kSatUpp_numpy)
    
    #Flatten Input Data:
    kSatUpp_numpy=kSatUpp_numpy.flatten()
    kSatLow_numpy=kSatLow_numpy.flatten()
    satVolMoistContUpp_numpy=satVolMoistContUpp_numpy.flatten()
    satVolMoistContLow_numpy=satVolMoistContLow_numpy.flatten()
    resVolMoistContUpp_numpy=resVolMoistContUpp_numpy.flatten()
    resVolMoistContLow_numpy=resVolMoistContLow_numpy.flatten()
    storCapUpp=storCapUpp.flatten()
    storUpp_numpy=storUpp_numpy.flatten()
    storLow_numpy=storLow_numpy.flatten()
    
    #Number of cells
    self.nCells=len(kSatUpp_numpy)
    
    #MODEL DIMENSIONS
    global dz, z, n, t, dt, timeunit
    #Depth:
    dzUpp_numpy=thickUpp.flatten()
    dzLow_numpy=thickLow.flatten()
    
    #ADDED FOR MULTIPLE LAYERS: START
    dz_pcr=np.vstack((dzLow_numpy, dzUpp_numpy)) #[m]; Thickness of soil layers in PCR-GlobWB.
    #Added for Layer Boundaries
    z_Tot=dz_pcr.sum(axis=0) #[m]; Total thickness of all soil layers together.
    
    #Next part could be useful if we want to divide the layers into a specific number.
    if self.layerFactorRichards>0:
        self.layerfractions=np.zeros((self.numberOfSoilLayers,int(self.layerFactorRichards*self.numberOfSoilLayers),self.nCells))
        for i in range(2):
            for j in range(self.layerFactorRichards):
                self.layerfractions[i, i*self.layerFactorRichards+j,:]=1/self.layerFactorRichards
        
        dz=np.array([np.dot(self.layerfractions[:,:,i].T, dz_pcr[:,i]) for i in range(self.nCells)]).T
        z=np.zeros(self.nCells)
        for i in range(len(dz[:,0])):
            z=np.vstack((z,dz[:i+1,:].sum(axis=0)))
    
    else:
        BoundUpp=np.argmax((self.layerBoundariesRichards>float(np.nanmin(dzUpp_numpy))) & (self.layerBoundariesRichards<float(np.nanmax(dzUpp_numpy)))) #Index were boundary of the upper layer is.
        BoundUpp_all=(self.layerBoundariesRichards>thickUpp.min()) & (self.layerBoundariesRichards<thickUpp.max()) 
        BoundUpp_all[BoundUpp]=0 #Other boundaries in range of upper soil thickness are deleted.
        self.layerBoundariesRichards=self.layerBoundariesRichards[BoundUpp_all==False]
        
        BoundLow=np.argmax((self.layerBoundariesRichards>z_Tot.min()) & (self.layerBoundariesRichards<z_Tot.max())) #Index were boundary of the lower layer is
        BoundLow_all=(self.layerBoundariesRichards>z_Tot.min()) & (self.layerBoundariesRichards<z_Tot.max())
        BoundLow_all[BoundLow]=0 #Other boundaries in range of lower soil thickness are deleted.
        self.layerBoundariesRichards=self.layerBoundariesRichards[BoundLow_all==False]
        
        #Give every cell the appropriate thicknesses.
        self.layerBoundariesRichards=np.array([ self.layerBoundariesRichards,]*self.nCells).T
        self.layerBoundariesRichards[BoundUpp,:]=dzUpp_numpy
        self.layerBoundariesRichards[BoundLow,:]=z_Tot
        
        z= self.layerBoundariesRichards[-1,:]- self.layerBoundariesRichards
        z=np.flip(z, axis=0) #[m]; Depth of the boundaries (z[-1] is surface, z[0] is bottom)
        dz=z[1:,:]-z[:-1,:] #[m]; Thickness of the layers.
        
        self.layerfractions=np.zeros((self.numberOfSoilLayers,len(dz),self.nCells)) #[fractions]; Used for conversion between the 2 PCR-GlobWB layers and the flexible layers.
        for j in range(len(dz)):
            if j<BoundLow-BoundUpp:
                i=0
            else:
                i=1
            self.layerfractions[i, j,:]=dz[j,:]/dz_pcr[i,:]
        self.layerFactorRichards=BoundLow-BoundUpp #Number of layers in which the lowest PCR-GlobWB layer is divided.
        #ADDED FOR MULTIPLE LAYERS: STOP
    
    #Shape of the arrays (axis=0 is soil layers, axis=1 is number of cells.)
    n=np.array([len(dz), self.nCells])
    
    #Time:
    dt=self.timestepRichards #[h]; timestep for the Richards model
    tsteps=2 #min 2 in dt hours. Timestep at which the result is evaluated.
    t = np.linspace(0,dt/24,tsteps) #Number of days per integration of Richards Model
    #dt=t.max()*24 #Number of hours per timestep
    timeunit=dt*3600 #[s]
    
    #SOIL CONDITIONS:
    global p, thetaS, thetaR, theta_soil
    
    def GuelphLoamDrying():
      '''
      Van Genuchten parameters. Have to be replaced by values from the model in a later version.
      '''
      pars={}
      #pars['thetaR']=0.218
      #pars['thetaS']=0.520
      pars['alpha']=1.15
      pars['n']=2.03
      pars['m']=1-1/pars['n']
      #pars['Ks']=0.316
      pars['neta']=0.5
      pars['Ss']=0.000001
      return pars
    
    #Give every cell the appropriate soil characteristics.
    pLow=np.array([]) #Characteristics of lower layer.
    for i in range(self.nCells):
        a=np.array([GuelphLoamDrying()])
        a[0]['Ks']=kSatLow_numpy[i]
        a[0]['thetaS']=satVolMoistContLow_numpy[i]
        a[0]['thetaR']=resVolMoistContLow_numpy[i]
        pLow=np.hstack((pLow,a))
    
    pUpp=np.array([]) #Characteristics of upper layer.
    for i in range(self.nCells):
        b=np.array([GuelphLoamDrying()])
        b[0]['Ks']=kSatUpp_numpy[i]
        b[0]['thetaS']=satVolMoistContUpp_numpy[i]
        b[0]['thetaR']=resVolMoistContUpp_numpy[i]
        pUpp=np.hstack((pUpp,b))
    p=np.vstack(([pLow,]*self.layerFactorRichards, [pUpp,]*(len(dz)-self.layerFactorRichards))) #Dictionary with soil characteristics.

        
    thetaS=np.vstack(([satVolMoistContLow_numpy,]*self.layerFactorRichards, [satVolMoistContUpp_numpy,]*(len(dz)-self.layerFactorRichards))) #[-]; Saturated moisture content of every layer and cell.
    theta_soil=1.-thetaS #[-]; Soil fraction
    
    thetaR=np.vstack(([resVolMoistContLow_numpy,]*self.layerFactorRichards, [resVolMoistContUpp_numpy,]*(len(dz)-self.layerFactorRichards))) #[-]; Residual moisture content of every layer and cell.
    
        
    #VEGETATION TYPES: 
    global veg, alpha #, Kcrop, ra_factor, evap_fraction, rootfraction
    def ReferenceGrass():
      pars={}
      pars['alpha']=0.25 #[-] Source:Wikipedia
      
      #Other characteristics could be used in future versions.
      #pars['rootdepth']=.3 #[m] 
      #pars['cc']=1 #[-] Crop coefficient: now only for midseason. Change this?
      #pars['ra_factor']=208 #From fao.org
      #pars['evap_fraction']=1 #Evaporation fraction: how much of the water is extracted from the soil?
      
      return pars
    
    veg=ReferenceGrass() #Set vegetation
    veg=np.repeat(veg, self.nCells)
    
    alpha=np.array([v['alpha'] for v in veg]) #[-]; Albedo
    
    #Other characteristics could be used in future versions.
    #Kcrop=np.array([v['cc'] for v in veg]) #Crop coefficients, based on Corbari et al., 2017 and FAO.
    #ra_factor=np.array([v['ra_factor'] for v in veg]) #Ra_factor, loosely based on FAO and Introduction to Physical Hydrology, Martin Hendriks.
    #evap_fraction=np.array([v['evap_fraction'] for v in veg]) #Evaporation fraction: how much of the evaporation is evaporated from the soil?
    
    #CONSTANTS
    global rho_a, cp, rho_s, cs, l_sat, l_dry, sigma, cw, dT, labda, rho_w
    rho_a=1.2 #[kg/m^3] Density of air.
    cp=1012 #[J/kg/K] Specific heat of air (isobaric)
    
    rho_s=1600 #[kg/m^3] Density for dry sand, Python Handbook
    cs=800 #[J/kg/K] Specific heat for dry sand, Python Handbook
    
    labda=2.5*1e6 #[J/kg]; Specific latent heat
    rho_w=1000 #[kg/m^3]; density of water
    
    l_sat=1.38*np.ones((len(z)-1, self.nCells)) #[W/m^2/K] for Silt Loam from Lu and Dong, 2015. --> Update this in a later version!
    l_dry=0.218*np.ones((len(z)-1, self.nCells))  #[W/m^2/K] for Silt Loam from Lu and Dong, 2015. --> Update this in a later version!
    
    sigma=5.67*1e-8 #Wm^-2K^-4; Stephan-Boltzmann constant
    cw=4186 #
    
    dT=0.2 #[K] Temperature step, for Newton-Raphson
    
    #INITIAL VALUES AND BOUNDARY CONDITIONS
    #IV TEMPERATURE
    #Initial soil temperature
    self.soiltemperature=280.*np.ones((n[0], self.nCells, 1)) #[K]; Array with soil temperature
    
    
    #IV RICHARDS
    global qBot, psiTop, psiBot
    qBot=np.array([None,]*self.nCells) #[m/day]; Fixed bottom flux (None for free drainage)
    psiTop=np.array([None,]*self.nCells) #[m]; Fixed matrix potential at the top. If there is water on the top; eg. 0.3.
    psiBot=np.array([None,]*self.nCells) #[m]; Fixed matrix potential at the bottom.
    
    
    storTot=np.vstack((storLow_numpy, storUpp_numpy)) #Initial storage in the subsurface
    self.storTotRich=np.array([np.dot(self.layerfractions[:,:,i].T, storTot[:,i]) for i in range(self.nCells)]).T #Coversion to multiple layers.
    
    
    #IV ENERGY BALANCE
    global dt_short, t_short, timeunit_short, SW_fractions, Tair_list
    dt_short=3 #[h]; time steps of energy balance
    t_short=np.linspace(0, dt_short/24, 2) #[day] Array for evaluation of energy balance
    timeunit_short=dt_short*3600 #[s] 
    SW_fractions=np.array([[0, 0, 0.15, 0.35, 0.35, 0.15, 0, 0]]).T #[%] Assumed daily radiation cycle.
    Tair_list=np.array([[-2.5, -2.5, 0, 0, 2.5, 2.5, 0, 0]]) #[K] Assumed daily air temperature cycle.
    
    msg='Initialisation Richards model finished!'
    logger.info(msg)
    return self
    

def runRichards(self, meteo, groundwater, currTimeStep):
    
    #CONVERSION TO NUMPY
    msg='Starting run with Richards model'
    logger.info(msg)
    
    Tair=pcr.pcr2numpy(meteo.temperature, np.nan)
    cloudCover=pcr.pcr2numpy(meteo.cloudCover, np.nan)
    vaporPressure=pcr.pcr2numpy(meteo.vaporPressure, np.nan)
    rsw=pcr.pcr2numpy(meteo.rsw, np.nan)
    
    storGroundwater=pcr.pcr2numpy(groundwater.storGroundwater, np.nan)
    
    infiltration_numpy=pcr.pcr2numpy(self.infiltration, 0.)
    storUpp_numpy=pcr.pcr2numpy(self.storUpp, np.nan)
    storLow_numpy=pcr.pcr2numpy(self.storLow, np.nan)
    actTranspiUpp_numpy=pcr.pcr2numpy(self.actTranspiUpp, np.nan)
    actTranspiLow_numpy=pcr.pcr2numpy(self.actTranspiLow, np.nan)
    actBareSoilEvap_numpy=pcr.pcr2numpy(self.actBareSoilEvap, np.nan)
    
    msg='Converting to numpy completed! We can start to run the model.'
    logger.info(msg)
    
    #%% INITIALISATION
    
    #Flatten the input
    infiltration_numpy=infiltration_numpy.flatten()
    actTranspiUpp_numpy=actTranspiUpp_numpy.flatten()
    actTranspiLow_numpy=actTranspiUpp_numpy.flatten()
    actBareSoilEvap_numpy=actBareSoilEvap_numpy.flatten()
    storUpp_numpy=storUpp_numpy.flatten()
    storLow_numpy=storLow_numpy.flatten()
    storGroundwater=storGroundwater.flatten()
    
    
    cloudCover=cloudCover.flatten()
    vaporPressure=vaporPressure.flatten()
    S0=rsw.flatten()
    Tair=Tair.flatten()
    Tair+=273.15 #To convert to Kelvin
    
    #Check NaN index:
    #Some cells are masked or do not take part in the model. These are not included in the computation.
    nanindexUpp=np.isnan(storUpp_numpy)
    nanindexLow=np.isnan(storLow_numpy)
    if (nanindexUpp==nanindexLow).all():
        msg='nanindex Upp == Low'
        logger.info(msg)
    else:
        msg='ERROR!!... nanindexUpp!=Low'
        logger.warning(msg)
        
    nanindex=np.nonzero(nanindexUpp)[0]
    nanindex=nanindex-np.arange(len(nanindex)) #The index changes if more values are inserted. This is the correction for this.
    
    #New n corrected for NaN-values.
    n1=n.copy()
    n1[1]=n1[1]-np.sum(nanindexUpp)
    
    
    #Remove NaN values:
    Tair=Tair[nanindexUpp==False]
    vaporPressure=vaporPressure[nanindexUpp==False]
    S0=S0[nanindexUpp==False]
    cloudCover=cloudCover[nanindexUpp==False]

    
    #COMPUTE MODEL INPUT:
    
    #Surface Energy Balance Input
    ra=69 #For Grass--> update this and make in mobile in a later version.
       
    RH=ComputeRH(vaporPressure, Tair)#[]; Relative humidity--> update this in a later version.
    RH[RH>1.]=0.99 #Since different data to compute the relative humidity, sometimes it reaches unrealistic values.
        
    S0=S0*(1-alpha[nanindexUpp==False])*24*3600 #[J/m^2]; Incoming short wave solar radiation in a day.
    
    plantevap=np.vstack((actTranspiLow_numpy[nanindexUpp==False], actTranspiUpp_numpy[nanindexUpp==False]+actBareSoilEvap_numpy[nanindexUpp==False]))
    plantevap=np.array([np.dot(self.layerfractions[:,:,i].T, plantevap[:,i]) for i in range(n1[1])]).T #[m]; Evaporation per layer
    #DOES THIS NOT TAKE THE WRONG LAYER FRACTIONS??
    ETpm=plantevap.sum(axis=0) #[m]; Total evaporation per cell
    
    #Added for 3H
    Ts=np.ones((len(SW_fractions), sum(nanindexUpp==False)))*Tair #[K]; Skin temperature
    Rn=np.zeros((len(SW_fractions), sum(nanindexUpp==False))) #[J/m^2]; Net Radiation
    H=np.zeros((len(SW_fractions), sum(nanindexUpp==False))) #[J/m^2]; Sensible Heat Flux
    G=np.zeros((len(SW_fractions), sum(nanindexUpp==False))) #[J/m^2]; Ground Heat Flux
    Sn=np.zeros((len(SW_fractions), sum(nanindexUpp==False))) #[J/m^2]; Net Incoming Shortwave Radiation
    Ln=np.zeros((len(SW_fractions), sum(nanindexUpp==False))) #[J/m^2]; Net Outgoing Longwave Radiation
    LE=np.zeros((len(SW_fractions), sum(nanindexUpp==False))) #[J/m^2]; Latent Heat Flux.
    
    #Add a daily cycle to climatological input.
    S0_3H=S0*SW_fractions 
    Tair_3H=Tair+Tair_list.T
    ETpm_3H=ETpm*SW_fractions
    
    #Richards Equation Input
    storTot_old=np.array([np.dot(self.layerfractions[:,:,i]>0, self.storTotRich[:,i]) for i in range(self.nCells)]).T
    storTot_old_real=np.vstack((storLow_numpy, storUpp_numpy))
    diff=storTot_old_real-storTot_old #Check if there is a difference the input from PCR-GlobWB and my input. Sometimes there is a small change--> is water extracted somewhere else?
    diff=np.array([np.dot(self.layerfractions[:,:,i].T, diff[:,i]) for i in range(self.nCells)]).T
    
    #Correct for possible different initial conditions.
    self.storTotRich+=diff #Difference is added, to correct for this.
    theta0=self.storTotRich[:,nanindexUpp==False]/dz[:,nanindexUpp==False] #[-]; Soil moisture content.
    theta0=theta0+thetaR[:, nanindexUpp==False]
    
    #Correct for possible oversaturated layers.
    excess=np.zeros(n1) #[m]; spill from oversaturated layers. Some layers were initially oversaturated due to the addition of diff.
    oversaturated=theta0>thetaS[:, nanindexUpp==False] #Index of oversaturated layers.
    excess[oversaturated]+=self.storTotRich[:, nanindexUpp==False][oversaturated]-(thetaS[:, nanindexUpp==False][oversaturated]-thetaR[:, nanindexUpp==False][oversaturated])*dz[:, nanindexUpp==False][oversaturated]
    self.storTotRich[:, nanindexUpp==False][oversaturated]=(thetaS[:, nanindexUpp==False][oversaturated]-thetaR[:, nanindexUpp==False][oversaturated])*dz[:, nanindexUpp==False][oversaturated]
    theta0[oversaturated]=thetaS[:, nanindexUpp==False][oversaturated]
    
    qTop=-infiltration_numpy #[m/day] Infiltration flux.
    
    #Define empty arrays
    interflow=np.zeros((int(n1[0]/2),n1[1])) #[m]
    satexcess=np.zeros((int(n1[0]/2),n1[1])) #[m]
        
    theta_run=theta0 #Soil moisture content for during the model, to keep the initial conditions for comparison.
        
    #Heat Function Input
    T0=self.soiltemperature[:,nanindexUpp==False, -1].reshape(n1) #[K]; Initial soil temperature
    
    
    #%% RUN MODEL
    '''
    SURFACE ENERGY BALANCE
    '''
    
    psi0=inv_thetaFun(theta0, p[:,nanindexUpp==False]) #[m]; Compute matric potential
    psi0=psi0.flatten() #Fluxmodel needs flattened input.
    q =fluxModel2(psi0,n1,p[:,nanindexUpp==False],qTop[nanindexUpp==False],qBot, psiTop[nanindexUpp==False],psiBot[nanindexUpp==False], z[:,nanindexUpp==False], T0); #Compute fluxes before model run.
    qOutUpp=q[self.layerFactorRichards,:] #[m/day]; Outflow upper layer at specific time steps (check for percUpp)
    qOutLow=q[0,:] #[m/day]; Outflow lower layer at specific time steps (check for percLow)
    psi0=psi0.reshape(n1)
    
    i_last=0 #Index needed for the second loop.
    
    for i in range(len(SW_fractions)): #Loop for Energy Balance
        Ts[i,:]=optimize.newton(func=radiation_f2, x0=Ts[i,:], fprime=radiation_deriv2, args=(dT, Tair_3H[i,:], T0[-1,:], RH, cloudCover, S0_3H[i,:], ra, rho_a, cp, dt_short, dz[-1, nanindexUpp==False], theta_run, thetaS[:, nanindexUpp==False], l_sat[:, nanindexUpp==False], l_dry[:, nanindexUpp==False], ETpm_3H[i,:]), maxiter=150)
        Rn[i,:], H[i,:], LE[i,:], G[i,:], Sn[i,:], Ln[i,:]=ComputeFluxes2(Ts[i,:], Tair_3H[i,:], T0[-1,:], RH, cloudCover, S0_3H[i,:], ra, rho_a, cp, dt_short, dz[-1, nanindexUpp==False], theta_run, thetaS[:, nanindexUpp==False], l_sat[:, nanindexUpp==False], l_dry[:, nanindexUpp==False], ETpm_3H[i,:])
        
# =============================================================================
#         '''
#         HEAT EQUATION: could also be placed in main loop (commented out)
#         '''
#         T0=T0.flatten() #Input data odeint has to be flattened.
#         T1 = odeint(heatFun3,T0,t,args=(q, z[:,nanindexUpp==False], cw, cs, theta_run, theta_soil[:, nanindexUpp==False], l_sat[:,nanindexUpp==False], l_dry[:,nanindexUpp==False], n1, dt), mxstep=1000)
#         T1=T1[-1,:].reshape(n1)
#         T0=T1
# =============================================================================
        
        if int((i+1)*dt_short)%dt==0: #Loop for Richards Model
            msg='Run model at t={}h'.format(int((i+1)*dt_short))
            logger.info(msg)
            
            #Added for 3H
            T0[-1, :]=Ts[:i,:].mean(axis=0) #Give the skin temperature as input for upper layer.
            
            #Subtract Evaporation.
            theta_top_new=(self.storTotRich[:, nanindexUpp==False]-plantevap*SW_fractions[i_last:i+1].sum())/dz[:,nanindexUpp==False]
            i_last=i+1
            theta_top_new=theta_top_new+thetaR[:, nanindexUpp==False] #[-]; New soil moisture content
            
            if self.numberOfCoresRichards>1:
            #Multicorerprocessing!
            #Send everything to self.
                global psi0_mp, t_mp, p_mp, qTop_mp,qBot_mp, psiTop_mp, psiBot_mp,z_mp,T0_mp, n2 #Global variables to use for multiprocessing.
                
                psi0_mp=psi0   
                t_mp=t
                p_mp=p[:,nanindexUpp==False]
                qTop_mp=qTop[nanindexUpp==False]
                qBot_mp=qBot[nanindexUpp==False]
                psiTop_mp=psiTop[nanindexUpp==False]
                psiBot_mp=psiBot[nanindexUpp==False]
                z_mp=z[:,nanindexUpp==False]
                T0_mp=T0
                n2=(n1[0],1) #Multiprocessing works per cell!
                
                #Actual computation
                pool = mp.Pool(self.numberOfCoresRichards)
                psi1 = pool.map(solve, range(n1[1]))
                pool.close()
                
                psi1=np.asarray(psi1)
                psi=psi1[:,1,:].T #Give the output the correct shape.
           
            else:
                #In case only one core is provided, run the model normally.
                psi0=psi0.flatten()
                psi1=odeint(RichardsModel_dz2,psi0,t,args=(n1, p[:,nanindexUpp==False], qTop[nanindexUpp==False],qBot,psiTop[nanindexUpp==False],psiBot[nanindexUpp==False],z[:,nanindexUpp==False], T0), mxstep=1000);
                psi=psi1[-1,:].reshape(n1)
                
            theta=thetaFun(psi, p[:,nanindexUpp==False]) #[-]; Compute new soil moisture content.
            
            '''
            HEAT EQUATION: can be placed in EB loop, but it makes it a little slower.
            '''
    
            T0=T0.flatten() #Input data odeint has to be flattened.
            T1 = odeint(heatFun3,T0,t,args=(q, z[:,nanindexUpp==False], cw, cs, theta, theta_soil[:, nanindexUpp==False], l_sat[:,nanindexUpp==False], l_dry[:,nanindexUpp==False], n1, dt), mxstep=1000)
            T1=T1[-1,:].reshape(n1)
            T0=T1
            
            
            theta=np.insert(theta, nanindex, np.nan, axis=1)
            
            #check for oversaturation:
            watercontent=theta[:,nanindexUpp==False]*dz[:,nanindexUpp==False]
            oversaturated=watercontent>(thetaS[:,nanindexUpp==False]*dz[:,nanindexUpp==False])
            excess[oversaturated]+=watercontent[oversaturated]-thetaS[:,nanindexUpp==False][oversaturated]*dz[:,nanindexUpp==False][oversaturated]
            theta[:,nanindexUpp==False][oversaturated]=thetaS[:,nanindexUpp==False][oversaturated]
            self.storTotRich=(theta-thetaR)*dz
            
            #Compute the new fluxes
            q =fluxModel2(psi,n1,p[:,nanindexUpp==False],qTop[nanindexUpp==False],qBot,psiTop[nanindexUpp==False],psiBot[nanindexUpp==False], z[:,nanindexUpp==False], T0); #[m/day]; Compute the new fluxes.
            #Make a list of the percolation (to check the waterbalance)
            qOutUpp=np.vstack((qOutUpp, q[self.layerFactorRichards,:]))
            qOutLow=np.vstack((qOutLow, q[0,:]))
            
            
            #New input for the next loop.
            theta_run=theta[:, nanindexUpp==False]
    
    #%% POST PROCESSING
    
    #Compute temperature difference Ts-Tair at 6am and 6pm (check daily cycle!)
    tempDeficit_6AM=Ts[1,:]-Tair_3H[1,:]
    tempDeficit_6PM=Ts[-3,:]-Tair_3H[-3,:]
    
    T1=np.insert(T1, nanindex, np.nan, axis=1)
     
    
    '''
    WB CHECK
    '''
    
    #ADDED FOR FIXING RICHARDS: START
    watercontent0=theta0*dz[:,nanindexUpp==False]-plantevap
    watercontent=theta[:,nanindexUpp==False]*dz[:,nanindexUpp==False]
    
    #Compute the percolation that must have taken place.
    percTot=np.zeros(n1)
    percTot[-1,:]=-(watercontent0[-1,:]-watercontent[-1,:]-qTop[nanindexUpp==False]-excess[-1,:])
    for layer in range(2, n1[0]+1):
        percTot[-layer,:]=-(watercontent0[-layer,:]-watercontent[-layer,:]-percTot[-layer+1,:]-excess[-layer,:])
    
    percUpp=percTot[self.layerFactorRichards,:]
    percLow=percTot[0,:]
    
    
# =============================================================================
#     #Weird Error: check if groundwater does not become empty..
#     no_gw=storGroundwater[nanindexUpp==False]-percLow<0
#     if no_gw.any():
#         msg='Groundwater is empty in some cells!'
#         logger.warning(msg)
#         diff=percLow[no_gw]-storGroundwater[nanindexUpp==False][no_gw]
#         percLow[no_gw]-=diff
#         j=0
#         while np.array(diff>0).any(): 
#             self.storTotRich[j,nanindexUpp==False][no_gw]-=diff
#             no_gw=self.storTotRich[j,nanindexUpp==False]<0
#             if np.array(no_gw).any():
#                 diff=-self.storTotRich[j,nanindexUpp==False][no_gw]
#                 self.storTotRich[j,nanindexUpp==False][no_gw]=0
#                 j+=1
#             else: 
#                 diff=0
# =============================================================================
            
    percUpp=np.insert(percUpp, nanindex, np.nan)
    percLow=np.insert(percLow, nanindex, np.nan)
    
    #Compute saturation excess and infiltration.
    interflow=excess[:self.layerFactorRichards,:]
    satexcess=excess[self.layerFactorRichards:,:]
    if self.layerFactorRichards>1:
        interflow=interflow.sum(axis=0)
    if (len(dz)-self.layerFactorRichards)>1:
        satexcess=satexcess.sum(axis=0)
    
    interflow=np.insert(interflow, nanindex, np.nan)
    satexcess=np.insert(satexcess, nanindex, np.nan)
    #ADDED FOR FIXING RICHARDS: STOP
    
    #Translate the new storage to the two layer PCR_GlobWB storage
    storTot_new=np.array([np.dot(self.layerfractions[:,:,i]>0, self.storTotRich[:,i]) for i in range(self.nCells)]).T
    storUpp_new=storTot_new[1,:]
    storLow_new=storTot_new[0,:]
    
    #Check again if there is no oversaturation...
    oversaturated=storUpp_new>storCapUpp
    if oversaturated.any():
        print('OVERSATURATED!')
    satexcess[oversaturated]+=storUpp_new[oversaturated]-storCapUpp[oversaturated]
    storUpp_new[oversaturated]=storCapUpp[oversaturated]
    
    #Translate the new temperature to the two layer PCR_GlobWB storage
    SoilTemp=np.array([np.dot(self.layerfractions[:,:,i], T1[:,i]) for i in range(self.nCells)]).T #Take the average over all layers in specific fractions for the temperature in both cells.
    SoilTempUpp=SoilTemp[1,:]
    SoilTempLow=SoilTemp[0,:]

     
    #Take the total daily energy flux.
    Rn=np.sum(Rn, axis=0)
    Sn=np.sum(Sn, axis=0)
    Ln=np.sum(Ln, axis=0)
    H=np.sum(H, axis=0)
    G=np.sum(G, axis=0)
    
    #Insert nan-values back in the right place.
    tempDeficit_6AM=np.insert(tempDeficit_6AM, nanindex, np.nan)
    tempDeficit_6PM=np.insert(tempDeficit_6PM, nanindex, np.nan)
    Rn=np.insert(Rn, nanindex, np.nan)
    Sn=np.insert(Sn, nanindex, np.nan)
    Ln=np.insert(Ln, nanindex, np.nan)
    H=np.insert(H, nanindex, np.nan)
    G=np.insert(G, nanindex, np.nan)
    ETpm=np.insert(ETpm, nanindex, np.nan)
    plantevap=np.insert(plantevap, nanindex, np.nan, axis=1)

    #Give everything the right shape.
    storUpp_new=storUpp_new.reshape(self.basin_shape)
    storLow_new=storLow_new.reshape(self.basin_shape)
    percUpp=percUpp.reshape(self.basin_shape)
    percLow=percLow.reshape(self.basin_shape)
    interflow=interflow.reshape(self.basin_shape)
    satexcess=satexcess.reshape(self.basin_shape)
    
    SoilTempUpp=SoilTempUpp.reshape(self.basin_shape)
    SoilTempLow=SoilTempLow.reshape(self.basin_shape)
    tempDeficit_6AM=tempDeficit_6AM.reshape(self.basin_shape)
    tempDeficit_6PM=tempDeficit_6PM.reshape(self.basin_shape)
    
    Rn=Rn.reshape(self.basin_shape)
    H=H.reshape(self.basin_shape)
    G=G.reshape(self.basin_shape)
    Ln=Ln.reshape(self.basin_shape)
    Sn=Sn.reshape(self.basin_shape)
    ETpm=ETpm.reshape(self.basin_shape)
    actBareSoilEvap_numpy=actBareSoilEvap_numpy.reshape(self.basin_shape)
    
    #Recompute the evaporation...
    plantevap=np.array([np.dot(self.layerfractions[:,:,i]>0, plantevap[:,i]) for i in range(self.nCells)]).T
    plantevapUpp=plantevap[1,:].reshape(self.basin_shape)-actBareSoilEvap_numpy
    plantevapLow=plantevap[0,:].reshape(self.basin_shape)
    
    #Translate to Energy
    LE=ETpm*labda*rho_w
    
    #Give everything the unit W/m^2
    timeunit1=24*3600
    Rn=Rn/timeunit1
    H=H/timeunit1
    G=G/timeunit1
    Ln=Ln/timeunit1
    LE=LE/timeunit1
    
    from scipy.stats import describe
    '''
    This part here was mainly for debugging. Could be used as future logging statements.
    '''
# =============================================================================
#     print('The New Ts:')
#     print(describe(Ts.mean(axis=0)))
#     
#     print('Evap difference:')
#     planttotal=np.sum(plantevap, axis=0)
#     print(describe(ETpm.flatten()[nanindexUpp==False]-planttotal[nanindexUpp==False], nan_policy='omit'))
#     
#     #ADDED FOR FIXING RICHARDS
#     print('Water Balance gap:')
#     WB=storUpp_numpy.flatten()[nanindexUpp==False]+storLow_numpy.flatten()[nanindexUpp==False]-storUpp_new.flatten()[nanindexUpp==False]-storLow_new.flatten()[nanindexUpp==False]-qTop[nanindexUpp==False]-planttotal[nanindexUpp==False]+percLow.flatten()[nanindexUpp==False]-interflow.flatten()[nanindexUpp==False]-satexcess.flatten()[nanindexUpp==False]
#     print(describe(WB,nan_policy='omit'))
#     
#     print('Water Balance Upp:')
#     WB_Upp=storUpp_numpy.flatten()[nanindexUpp==False]-storUpp_new.flatten()[nanindexUpp==False]-qTop[nanindexUpp==False]-plantevapUpp.flatten()[nanindexUpp==False]-actBareSoilEvap_numpy.flatten()[nanindexUpp==False]+percUpp.flatten()[nanindexUpp==False]-satexcess.flatten()[nanindexUpp==False]
#     print(describe(WB_Upp,nan_policy='omit'))
#     
#     print('Water Balance Low:')
#     WB_Low=storLow_numpy.flatten()[nanindexUpp==False]-storLow_new.flatten()[nanindexUpp==False]-percUpp.flatten()[nanindexUpp==False]-plantevapLow.flatten()[nanindexUpp==False]+percLow.flatten()[nanindexUpp==False]-interflow.flatten()[nanindexUpp==False]
#     print(describe(WB_Low,nan_policy='omit'))
#     
#     print('StorCapUpp: exceeded?')
#     SC_Upp=storUpp_new.flatten()[nanindexUpp==False]-storCapUpp[nanindexUpp==False]
#     print(describe(SC_Upp,nan_policy='omit'))
# =============================================================================
    
    #Check if the found percolation is physical.
    qOutUpp=qOutUpp[1:,:]
    qOutLow=qOutLow[1:,:]
    tol=1e-6
    indexUpp=(percTot[self.layerFactorRichards,:]<=qOutUpp.max(axis=0)+tol)*(qOutUpp.min(axis=0)-tol<=percTot[self.layerFactorRichards,:])
    if (indexUpp==False).any():
        msg='WB ERROR IN PERCOLATION UPP!'
        logger.debug(msg)
        
    indexLow=(qOutLow.min(axis=0)-tol<=percTot[0,:]) * (percTot[0,:]<=qOutLow.max(axis=0)+tol)
    if (indexLow==False).any():
        msg='WB ERROR IN PERCOLATION LOW!'
        logger.debug(msg)
    
    msg='PercUpp Correct:{}; PercLow Correct:{}'.format(np.sum(indexUpp), np.sum(indexLow))
    logger.debug(msg)
    
    percLow=-percLow
    percUpp=-percUpp
                      
    #%% CONVERSION TO PCRASTER
    
    #Update the new values:
    self.storUpp=pcr.numpy2pcr(pcr.Scalar, np.asarray(storUpp_new.tolist()), np.nan) #To list conversion is needed since array is masked or something? Has to be fixed in next version!
    self.storLow=pcr.numpy2pcr(pcr.Scalar, np.asarray(storLow_new.tolist()), np.nan)
    self.percUpp=pcr.numpy2pcr(pcr.Scalar, np.asarray(percUpp.tolist()), np.nan)
    self.percLow=pcr.numpy2pcr(pcr.Scalar, np.asarray(percLow.tolist()), np.nan)

    self.soilTempUpp=pcr.numpy2pcr(pcr.Scalar, SoilTempUpp, np.nan)
    self.soilTempLow=pcr.numpy2pcr(pcr.Scalar, SoilTempLow, np.nan)
    self.tempDeficit_6AM=pcr.numpy2pcr(pcr.Scalar, tempDeficit_6AM, np.nan)
    self.tempDeficit_6PM=pcr.numpy2pcr(pcr.Scalar, tempDeficit_6PM, np.nan)
    
    self.longWaveRad=pcr.numpy2pcr(pcr.Scalar, Ln, np.nan)
    self.netSW=pcr.numpy2pcr(pcr.Scalar, Sn, np.nan)
    self.netRad=pcr.numpy2pcr(pcr.Scalar, Rn, np.nan)
    self.latentHF=pcr.numpy2pcr(pcr.Scalar, LE, np.nan)
    self.sensibleHF=pcr.numpy2pcr(pcr.Scalar, H, np.nan)
    self.groundHF=pcr.numpy2pcr(pcr.Scalar, G, np.nan)
    self.soiltemperature=np.expand_dims(T1, axis=2)
    self.actTranspiUpp=pcr.numpy2pcr(pcr.Scalar, np.asarray(plantevapUpp.tolist()), np.nan)
    self.actTranspiLow=pcr.numpy2pcr(pcr.Scalar, np.asarray(plantevapLow.tolist()), np.nan)
    self.actBareSoilEvap=pcr.numpy2pcr(pcr.Scalar, np.asarray(actBareSoilEvap_numpy.tolist()), np.nan)
    
    self.interflow=pcr.numpy2pcr(pcr.Scalar, np.asarray(interflow.tolist()), np.nan)
    self.satExcess=pcr.numpy2pcr(pcr.Scalar, np.asarray(satexcess.tolist()), np.nan)
     
    #Set everything to zero what has not yet been included in the model.
    self.capRiseUpp=pcr.min(0., self.capRiseUpp)
    
    #Keep positive and negative fluxes seperate:
    self.capRiseLow=pcr.min(0., self.percLow)*-1.
    self.percLow=pcr.max(0., self.percLow)
    
    msg='Succesfull Run!'
    logger.info(msg)    
                
    return self

    
    
    
#ADDED FOR MULTIPROCESSING: START
def solve(i):
    '''
    Works only for 2 layers system
    '''
    #n2=(2,1)
    psi2=psi0_mp[:,i].flatten()
    z_solve=np.expand_dims(z_mp[:,i],axis=1)
    T0_solve=np.expand_dims(T0_mp[:,i],axis=1)
    p_solve=np.expand_dims(p_mp[:,i],axis=1)
    
    sol = odeint(RichardsModel_dz2,psi2,t_mp,args=(n2, p_solve, qTop_mp[i],qBot_mp[i], psiTop_mp[i],psiBot_mp[i],z_solve, T0_solve), mxstep=3500);
    return sol
#ADDED FOR MULTIPROCESSING: STOP                


def RichardsModel_dz2(psi,t,n,p,qTop,qBot,psiTop,psiBot, z, T):
 
    psi = psi.reshape(n)
    zMid = np.hstack(z[:-1, :] + (z[1::,:]-z[:-1,:])/2).reshape(n) #Compute midpoint of each cell
    dz = np.abs(np.hstack([zMid[0,:]-z[0,:], (zMid[1::,:]-zMid[:-1,:]).flatten(), z[-1,:]-zMid[-1,:]])).reshape((n[0]+1,n[1])) #Compute distance between midpoint (location of fluxes)
    
    
    # Basic properties:
    C=CFun(psi,p)

    # initialize vectors:
    q=np.zeros((n[0]+1, n[1]))

    # Upper boundary
    if np.isnan(qTop).any():
        flooded=np.isnan(qTop)
        KTop=KFun1(np.zeros(np.sum(flooded))+psiTop[flooded],p[n[0]-1, flooded], T[n[0]-1, flooded])
        q[n[0],flooded]=-KTop*((psiTop[flooded]-psi[n[0]-1,flooded])/(dz[-1,flooded]*2)+1)
        q[n[0],np.isfinite(qTop)]=qTop[np.isfinite(qTop)]
        #print(q[n[0], :])
    else:
        q[n[0],:]=qTop

    # Lower boundary
    if np.array([qBot==None]).all():
        #if psiBot is None:
        if np.array([psiBot==None]).all():
            # Free drainage
            KBot=KFun1(np.zeros(n[1])+psi[0,:],p[0], T[0])
            q[0,:]=-KBot
        else:
            # Type 1 boundary
            KBot=KFun1(np.zeros(1)+psiBot,p[0], T[0])
            q[0,:]=-KBot*((psi[0]-psiBot)/(dz[0,:]*2)+1.0)
    else:
        # Type 2 boundary
        q[0,]=qBot

    # Internal nodes
    i=np.arange(0,n[0]-1)
    Knodes=KFun1(psi,p, T)
    Kmid=(Knodes[i+1,:]+Knodes[i,:])/2.0
    j=np.arange(1,n[0])
    
    q[j]=-Kmid*((psi[i+1,:]-psi[i,:])/dz[j,:]+1.0)

    # Continuity
    i=np.arange(0,n[0])
    
    #ADDED FOR PERCOLATION
    dz1=z[1:,:]-z[:-1,:] #For psi the thickness of entire layer is needed!
    dpsidt=((-(q[i+1,:]-q[i,:])/dz1[i,:])/C).flatten().astype(float)
    
    return dpsidt


def heatFun3(T, t, q, z, Cw, Cs, thetaW, thetaS, l_sat, l_dry, n, dt):
    T = T.reshape(n)
    l=ThermalConduct(thetaW[1:-1], thetaS[1:-1], l_sat[1:-1], l_dry[1:-1], T[1:-1])
    l*=dt*3600
    heat_cap=heatcapacity1(thetaW[1:-1], thetaS[1:-1], T[1:-1])

    zMid = np.hstack(z[:-1, :] + (z[1::,:]-z[:-1,:])/2).reshape(n) #Compute midpoint of each cell
    dz = np.abs(np.hstack([zMid[0,:]-z[0,:], (zMid[1::,:]-zMid[:-1,:]).flatten(), z[-1,:]-zMid[-1,:]])).reshape((n[0]+1,n[1])) #Compute distance between midpoint (location of fluxes)
    
    
    Ew = (T)*Cw*thetaW
    Es = (T)*Cs*thetaS
    E = Ew + Es
    i=np.arange(1,n[0]-1)
    
    deltaT = T[i+1]-2*T[i]+T[i-1]
    deltaE = 0.5*(Ew[i-1]-Ew[i])/dz[i-1] + 0.5*(Ew[i]-Ew[i+1])/dz[i]
    
    conduc = deltaT/(dz[i-1]*dz[i]) *l/heat_cap #Conduction
    convec = deltaE * -q[i]/heat_cap #Convection
    dTdt = np.vstack([conduc-convec,np.zeros(n[1])]).flatten()
    dTdt = np.hstack([dTdt[:n[1]], dTdt]) #Add this give the lower layer the same derivative as the layer above it!
    return dTdt


def fluxModel2(psi,n,p,qTop,qBot,psiTop,psiBot, z, T):
    psi = psi.reshape(n)
    
    zMid = np.hstack(z[:-1, :] + (z[1::,:]-z[:-1,:])/2).reshape(n) #Midpoint of each layer
    dz = np.abs(np.hstack([zMid[0,:]-z[0,:], (zMid[1::,:]-zMid[:-1,:]).flatten(), z[-1,:]-zMid[-1,:]])).reshape((n[0]+1,n[1])) #Distance between midpoints.
    
    # initialize vectors:
    q=np.zeros((n[0]+1, n[1]))

    # Upper boundary
    if np.isnan(qTop).any():
        flooded=np.isnan(qTop)
        KTop=KFun1(np.zeros(np.sum(flooded))+psiTop[flooded],p[n[0]-1, flooded], T[n[0]-1, flooded])
        q[n[0],flooded]=-KTop*((psiTop[flooded]-psi[n[0]-1,flooded])/(dz[-1,flooded]*2)+1)
        q[n[0],np.isfinite(qTop)]=qTop[np.isfinite(qTop)]
    else:
        q[n[0],:]=qTop

    # Lower boundary
    if np.array([qBot==None]).all():
        if np.array([psiBot==None]).all():
            # Free drainage
            KBot=KFun1(np.zeros(n[1])+psi[0,:],p[0], T[0])
            q[0,:]=-KBot
        else:
            # Type 1 boundary
            KBot=KFun1(np.zeros(1)+psiBot,p[0], T[0])
            q[0,:]=-KBot*((psi[0]-psiBot)/(dz[0,:]*2)+1.0)
    else:
        # Type 2 boundary
        q[0,]=qBot
        
    # Internal nodes
    i=np.arange(0,n[0]-1)
    Knodes=KFun1(psi,p,T)
    Kmid=(Knodes[i+1,:]+Knodes[i,:])/2.0
    j=np.arange(1,n[0])
    q[j]=-Kmid*((psi[i+1,:]-psi[i,:])/dz[j,:]+1.0)
    return q


#CHANGED FOR FIXING RICHARDS
def thetaFun(psi,pars):
    '''
    Formula 3 from Van Genuchten 1980.
    Returns theta
    '''
    constant=0
    if psi>=0.:
        Se=1.
        constant=CFun(psi, pars)*psi #To include oversaturation.
    else:
        Se=(1+abs(psi*pars['alpha'])**pars['n'])**(-pars['m'])
    return pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se+constant
 
thetaFun=np.vectorize(thetaFun, otypes=[float])


def inv_thetaFun(theta,pars):
    '''
    Comput psi from theta. Inverse thetaFun
    '''
    if np.isnan(theta):
        return np.nan
    if theta<=pars['thetaR']:
        theta=pars['thetaR']+1e-8 #To make sure it does not become imaginary
    Se=(theta-pars['thetaR'])/(pars['thetaS']-pars['thetaR'])
    if Se>=1:
        psi=0
    else:
        psi=-1*((Se**(-1/pars['m'])-1)**(1/pars['n'])/pars['alpha'])
    return psi

inv_thetaFun=np.vectorize(inv_thetaFun, otypes=[float]) #otypes=[np.ndarray])
    

def CFun(psi,pars):
    '''
    How much psi change as theta changes. Van Genuchten.
    '''
    if psi>=0.:
        Se=1.
    else:
        Se=(1+abs(psi*pars['alpha'])**pars['n'])**(-pars['m'])
    dSedh=pars['alpha']*pars['m']/(1-pars['m'])*Se**(1/pars['m'])*(1-Se**(1/pars['m']))**pars['m']
    return Se*pars['Ss']+(pars['thetaS']-pars['thetaR'])*dSedh

CFun = np.vectorize(CFun, otypes=[float])

def KFun1(psi,pars,T):
    '''
    Unsaturated Hydraulic Conductivity.
    '''
    if T<=273.15:
        return 0 #No flow if it freezes!
    if psi>=0.:
        Se=1.
    else:
        Se=(1+abs(psi*pars['alpha'])**pars['n'])**(-pars['m'])
    return pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2
  
KFun1 = np.vectorize(KFun1, otypes=[float])

def heatcapacity1(theta,thetaS, T):
    rho_a=1.2 #air density kg/m^3
    cp=1012 #[J/kg/K] Specific heat of air 
    rho_s=1600 #[kg/m^3] density for dry sand, Python Handbook
    cs=800 #[J/kg/K] specific heat for dry sand, Python Handbook
    if T<=273.15: #Include freezing!
        ci=2027 #[J/kg/K]Specifici heat ice; source: https://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html
        rho_i=927#[kg/m^2] density ice
        return ci*rho_i*theta+cp*rho_a*(thetaS-theta)+cs*rho_s*(1-thetaS)
    cw=4186 #[J/kg/K] Specific heat water
    rho_w=1000 #[kg/m^-3]
    return cw*rho_w*theta+cp*rho_a*(thetaS-theta)+cs*rho_s*(1-thetaS)

heatcapacity1=np.vectorize(heatcapacity1, otypes=[float])

def ThermalConduct(theta_top, thetaS, l_sat, l_dry, T):
    '''
    Based on Jules model, simple linear.
    '''
    if T<=273.15: #Include freezing!
        l_ice=2.25 #[W/m/K]
        return (l_ice-l_dry)*(theta_top/thetaS)+l_dry
    Ks=(l_sat-l_dry)*(theta_top/thetaS)+l_dry
    #print(Ks)
    return Ks

ThermalConduct=np.vectorize(ThermalConduct, otypes=[float])


#EB FUNCTIONS
def Gflux(Ts, T_soil, timeunit, dz, Ks):
    '''
    Compute Ground Heat flux.
    '''
    G=Ks*(Ts-T_soil)/dz*timeunit
    return G


def longWave_out2(Ts, Tair, RH, cloudfraction, timeunit):
    '''
    Physical Hydrology, Dingman; p. 261.
    '''
    
    e_s=0.6108*np.exp(17.27*(Tair-273.15)/(237.3+(Tair-273.15))) 
    e_a=RH*e_s
    epsilon_a=(1-0.84*cloudfraction)*(0.83-0.18*np.exp(-1.54*e_a))+0.84*cloudfraction #Emissivity air.
    
    epsilon_s=0.975 #Emissivity soil; https://www.engineeringtoolbox.com/radiation-heat-emissivity-d_432.html
    
    sigma=5.67*1e-8 #Wm^-2K^-4
    Ln=epsilon_s*epsilon_a*sigma*Tair**4-epsilon_s*sigma*Ts**4
    Ln=Ln*timeunit
    
    return -Ln


def computeH(Tair, Ts, ra, rho_a, cp, timeunit):
    rho_a=1.2 #air density kg/m^3
    cp=1012 #[J/kg/K] Specific heat of air 
    H=-rho_a*cp*(Tair-Ts)/(ra)*timeunit
    return H


def ComputeRH(vaporPressure,Tair):
    '''
    Teten's equation.
    '''
    e_s=0.6108*np.exp(17.27*(Tair-273.15)/(237.3+(Tair-273.15)))
    return vaporPressure/e_s

def ComputeFluxes2(Ts, T_air, T_soil, RH, cloudfraction, Sin, ra, rho_a, cp, dt, dz, theta, thetaS, l_sat, l_dry, LE):
    '''
    Returns Energy Fluxes in Energy Balance.
    '''
    #timeunit=dt*3600 #[s]
    Lout = longWave_out2(Ts, T_air, RH, cloudfraction, timeunit)

    H = computeH(T_air, Ts, ra, rho_a, cp, timeunit)
    Ks=ThermalConduct(theta[-1], thetaS[-1], l_sat[-1], l_dry[-1], T_soil)
    G=Gflux(Ts, T_soil, timeunit, dz, Ks)
    Rn = Sin - Lout
    #LE=LE*dt/24 #Give it the right unit
    
    return Rn, H, LE, G, Sin, Lout

def radiation_f2(Ts, dT, T_air, T_soil, RH, cloudfraction, Sin, ra, rho_a, cp, dt, dz, theta, thetaS, l_sat, l_dry, LE):
    '''
    Energy Balance (for Newton-Raphson)
    '''
    labda=2.5*1e6
    rho_w=1000 #kg/m^3
    Rn, H, LE, G, Sin, Lout=ComputeFluxes2(Ts, T_air, T_soil, RH, cloudfraction, Sin, ra, rho_a, cp, dt, dz, theta, thetaS, l_sat, l_dry, LE)
    
    return Rn-H-G-rho_w*labda*LE

def radiation_deriv2(Ts, dT, T_air, T_soil, RH, cloudfraction, Sin, ra, rho_a, cp, dt, dz, theta, thetaS, l_sat, l_dry, LE):
    '''
    Simple erivative Energy Balance (for Newton-Raphson)
    '''
    f=radiation_f2(Ts-dT, dT, T_air, T_soil, RH, cloudfraction, Sin, ra, rho_a, cp, dt, dz, theta, thetaS, l_sat, l_dry, LE)
    fdT=radiation_f2(Ts+dT,dT, T_air, T_soil, RH, cloudfraction, Sin, ra, rho_a, cp, dt, dz, theta, thetaS, l_sat, l_dry, LE)

    return (fdT-f)/(2*dT)


#FOR FUTURE USE: A SIMPLE SNOWCOVER MODEL: A FIRST START. IMPLEMENTED IN FullModel_Raam.py
    
# =============================================================================
# #Snow Cover
#     if Snow_Cover.any():
#         snowvariables=nb.snow_cover(Snow_Cover, snowdepth[count,:], T_soil, qTop1, veg, snowcount[count,:], alpha, dt_short)
#         for i in range(nCells):
#             snowdepth[count+1,i]=snowvariables[i][0]
#             snowcount[count+1,i]=snowvariables[i][1]
#             alpha[i]=snowvariables[i][2]
#             qTop2[i]=snowvariables[i][3]
#         alpha_list[count+1, :]=alpha
#         qTopmelt_list[count+1, :]=qTop2
#         qTop+=qTop2
#     else:
#         qTop+=qTop1  
# =============================================================================

# =============================================================================
# def snow_cover(cover, snowdepth, Tair, qTop, veg, snowcount, alpha, dt):
#     #snow count in timesteps
#     if cover==True:
#         alpha_fresh=0.9
#         alpha_old=veg['alpha']
#         lifetime=2*24/dt#dagen
#         if Tair<272 and qTop<0:
#                 alpha=alpha_fresh
#                 snowdepth1=-qTop+snowdepth
#                 qTop=0
#                 snowcount=0
#                 
#         elif Tair<277:
#             if snowdepth<0.0005:
#                 snowdepth1=0
#                 qTop-=snowdepth
#                 snowcount=0
#                 alpha=veg['alpha']
#             else:
#                 snowdepth1=snowdepth*np.exp(-snowcount/lifetime)
#                 qTop-=snowdepth-snowdepth1
#                 alpha=alpha_old+(alpha_fresh-alpha_old)*np.exp(-snowcount/lifetime)
#                 snowcount+=1
#         elif snowdepth:
#             snowdepth1=0
#             qTop-=snowdepth
#             snowcount=0
#             alpha=veg['alpha']
#         else:
#             snowdepth1=0
#             alpha=veg['alpha']
#     else:
#         snowdepth1=0
#         alpha=veg['alpha']
#     return snowdepth1, snowcount, alpha, qTop
# 
# snow_cover=np.vectorize(snow_cover, otypes=[np.ndarray], cache=False)
# 
# 
# def snow_temperature(snowdepth, Ts):
#     #Ts computation based on Rankinen et al., 2004 (Adaptation of the Integrated..)
#     CF=-0.025*100 #[m^-1] Rankinen et al., 2004
#     Tnew=(Ts-273.15)*np.exp(CF*snowdepth)+273.15
#     return Tnew
# 
# snow_temperature=np.vectorize(snow_temperature, otypes=[float], cache=False)
# =============================================================================
