import numpy as np
from spacetime_minkowski import SpaceTime_Minkowski

class SpaceTime_Schwarschild(SpaceTime_Minkowski):
   
   G = 1   # Newton's constant
   M = 1   # mass
   r_s = 2*G*M/SpaceTime_Minkowski.C**2

   def MetricTensor(self,r,theta,phi):

      r_s = self.r_s
      C = self.C

      g__dt_dt         = -(1-r_s/r)*C**2 if not r == 0 else np.inf
      g__dt_dr         = 0
      g__dt_dtheta     = 0
      g__dt_dphi       = 0
      g__dr_dr         = 1/(1-r_s/r) if not r == 0 else -np.inf
      g__dt_dtheta     = 0
      g__dr_dphi       = 0
      g__dtheta_dtheta = r**2
      g__dtheta_dphi   = 0
      g__dphi_dphi     = r**2*np.sin(theta)**2
      
      return g__dt_dt,g__dt_dr,g__dt_dtheta,g__dt_dphi,g__dr_dr, \
             g__dt_dtheta,g__dr_dphi,g__dtheta_dtheta,g__dtheta_dphi, \
             g__dphi_dphi

   def ChristoffelSymbols(self,r,theta,phi):

      r_s = self.r_s
      C = self.C

      Gamma_t__t_t         = 0
      Gamma_t__t_r         = r_s/(2*r*(r-r_s)) if not r == 0 else np.inf
      Gamma_t__t_theta     = 0
      Gamma_t__t_phi       = 0
      Gamma_t__r_r         = 0
      Gamma_t__r_theta     = 0
      Gamma_t__r_phi       = 0
      Gamma_t__theta_theta = 0
      Gamma_t__theta_phi   = 0
      Gamma_t__phi_phi     = 0

      Gamma_r__t_t         = C**2*r_s*(r-r_s)/(2*r**3) if not r == 0 else np.inf
      Gamma_r__t_r         = 0
      Gamma_r__t_theta     = 0
      Gamma_r__t_phi       = 0
      Gamma_r__r_r         = -r_s/(2*r*(r-r_s)) if not r == 0 else -np.inf
      Gamma_r__r_theta     = 0
      Gamma_r__r_phi       = 0
      Gamma_r__theta_theta = -(r-r_s)
      Gamma_r__theta_phi   = 0
      Gamma_r__phi_phi     = -(r-r_s)*np.sin(theta)**2

      Gamma_theta__t_t         = 0
      Gamma_theta__t_r         = 0
      Gamma_theta__t_theta     = 0
      Gamma_theta__t_phi       = 0
      Gamma_theta__t_phi       = 0
      Gamma_theta__r_r         = 0
      Gamma_theta__r_theta     = 1/r if not r == 0 else np.inf
      Gamma_theta__r_phi       = 0
      Gamma_theta__theta_theta = 0
      Gamma_theta__theta_phi   = 0
      Gamma_theta__phi_phi     = -np.sin(theta)*np.cos(theta)

      Gamma_phi__t_t         = 0
      Gamma_phi__t_r         = 0
      Gamma_phi__t_theta     = 0
      Gamma_phi__t_phi       = 0
      Gamma_phi__r_r         = 0
      Gamma_phi__r_theta     = 0
      Gamma_phi__r_phi       = 1/r if not r == 0 else np.inf
      Gamma_phi__theta_theta = 0
      Gamma_phi__theta_phi   = np.cos(theta)/np.sin(theta) if not theta==0 else np.inf
      Gamma_phi__phi_phi     = 0

      Gamma_t = np.array([Gamma_t__t_t     ,Gamma_t__t_r     ,Gamma_t__t_theta     ,Gamma_t__t_phi, \
                          Gamma_t__t_r     ,Gamma_t__r_r     ,Gamma_t__r_theta     ,Gamma_t__r_phi, \
                          Gamma_t__t_theta ,Gamma_t__r_theta ,Gamma_t__theta_theta ,Gamma_t__theta_phi, \
                          Gamma_t__t_phi   ,Gamma_t__r_phi   ,Gamma_t__theta_phi   ,Gamma_t__phi_phi])

      Gamma_r = np.array([Gamma_r__t_t     ,Gamma_r__t_r     ,Gamma_r__t_theta     ,Gamma_r__t_phi, \
                          Gamma_r__t_r     ,Gamma_r__r_r     ,Gamma_r__r_theta     ,Gamma_r__r_phi, \
                          Gamma_r__t_theta ,Gamma_r__r_theta ,Gamma_r__theta_theta ,Gamma_r__theta_phi, \
                          Gamma_r__t_phi   ,Gamma_r__r_phi   ,Gamma_r__theta_phi   ,Gamma_r__phi_phi])

      Gamma_theta = np.array([Gamma_theta__t_t     ,Gamma_theta__t_r     ,Gamma_theta__t_theta     ,Gamma_theta__t_phi, \
                              Gamma_theta__t_r     ,Gamma_theta__r_r     ,Gamma_theta__r_theta     ,Gamma_theta__r_phi, \
                              Gamma_theta__t_theta ,Gamma_theta__r_theta ,Gamma_theta__theta_theta ,Gamma_theta__theta_phi, \
                              Gamma_theta__t_phi   ,Gamma_theta__r_phi   ,Gamma_theta__theta_phi   ,Gamma_theta__phi_phi])

      Gamma_phi = np.array([Gamma_phi__t_t     ,Gamma_phi__t_r     ,Gamma_phi__t_theta     ,Gamma_phi__t_phi, \
                            Gamma_phi__t_r     ,Gamma_phi__r_r     ,Gamma_phi__r_theta     ,Gamma_phi__r_phi, \
                            Gamma_phi__t_theta ,Gamma_phi__r_theta ,Gamma_phi__theta_theta ,Gamma_phi__theta_phi, \
                            Gamma_phi__t_phi   ,Gamma_phi__r_phi   ,Gamma_phi__theta_phi   ,Gamma_phi__phi_phi])

      return Gamma_t,Gamma_r,Gamma_theta,Gamma_phi
