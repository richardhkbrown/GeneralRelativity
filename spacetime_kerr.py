import numpy as np
from spacetime_schwarschild import SpaceTime_Schwarschild

class SpaceTime_Kerr(SpaceTime_Schwarschild):

   J = 1 # angular momentum
   a = J/(SpaceTime_Schwarschild.M*SpaceTime_Schwarschild.C)

   def MetricTensor(self,r,theta,phi):

      r_s = self.r_s
      M = self.M
      a = self.a
      C = self.C

      Sigma = r**2 + a**2*np.cos(theta)**2
      Delta = r**2-r_s*r+a**2
      if Delta<0:
         raise TypeError('Delta %F has invalid value.'%Delta)      

      g__dt_dt         = -(1-r_s*r/Sigma)*C**2 if not Sigma == 0 else np.inf
      g__dt_dr         = 0
      g__dt_dtheta     = 0
      g__dt_dphi       = -2*r_s*a*r*np.sin(theta)**2/Sigma*C if not Sigma == 0 else np.inf
      g__dr_dr         = Sigma / Delta if not Delta == 0 else -np.inf
      g__dt_dtheta     = 0
      g__dr_dphi       = 0
      g__dtheta_dtheta = Sigma
      g__dtheta_dphi   = 0
      g__dphi_dphi     = (r**2+a**2+r_s*a**2*r*np.sin(theta)**2/Sigma)*np.sin(theta)**2 \
                          if not Sigma == 0 else np.inf
      
      return g__dt_dt,g__dt_dr,g__dt_dtheta,g__dt_dphi,g__dr_dr, \
             g__dt_dtheta,g__dr_dphi,g__dtheta_dtheta,g__dtheta_dphi, \
             g__dphi_dphi

   def ChristoffelSymbols(self,r,theta,phi):

      r_s = self.r_s
      a = self.a
      C = self.C

      Sigma = r**2 + a**2*np.cos(theta)**2
      Delta = r**2 - r_s*r + a**2
      A = (r**2+a**2)**2-a**2*Delta*np.sin(theta)**2

      Gamma_t__t_t         = 0
      Gamma_t__t_r         = r_s*(r**2+a**2)*(r**2-a**2*np.cos(theta)**2)/(2*Sigma**2*Delta) \
                             if not Sigma == 0 or Delta == 0 else np.inf
      Gamma_t__t_theta     = -r_s*a**2*r*np.sin(theta)*np.cos(theta)/Sigma**2 \
                             if not Sigma == 0 else -np.inf
      Gamma_t__t_phi       = 0
      Gamma_t__r_r         = 0
      Gamma_t__r_theta     = 0
      Gamma_t__r_phi       = r_s*a*np.sin(theta)**2*(a**2*np.cos(theta)**2*(a**2-r**2)-r**2*(a**2+3*r**2))/ \
                             (2*C*Sigma**2*Delta) if not Sigma == 0 or Delta == 0 else np.inf
      Gamma_t__theta_theta = 0
      Gamma_t__theta_phi   = r_s*a**3*r*np.sin(theta)**3*np.cos(theta)/(C*Sigma**2) \
                             if not Sigma == 0 else -np.inf
      Gamma_t__phi_phi     = 0

      Gamma_r__t_t         = C**2*r_s*Delta*(r**2-a**2*np.cos(theta)**2)/(2*Sigma**3) \
                             if not Sigma == 0 else -np.inf
      Gamma_r__t_r         = 0
      Gamma_r__t_theta     = 0
      Gamma_r__t_phi       = -C*Delta*r_s*a*np.sin(theta)**2*(r**2-a**2*np.cos(theta)**2)/(2*Sigma**3) \
                             if not Sigma == 0 else -np.inf
      Gamma_r__r_r         = (2*r*a**2*np.sin(theta)**2-r_s*(r**2-a**2*np.cos(theta)**2))/(2*Sigma*Delta) \
                             if not Sigma == 0 or Delta == 0 else np.inf
      Gamma_r__r_theta     = -a**2*np.sin(theta)*np.cos(theta)/Sigma \
                             if not Sigma == 0 else -np.inf
      Gamma_r__r_phi       = 0
      Gamma_r__theta_theta = -r*Delta/Sigma if not Sigma == 0 else np.inf
      Gamma_r__theta_phi   = 0
      Gamma_r__phi_phi     = Delta*np.sin(theta)**2/(2*Sigma**3)*(-2*r*Sigma**2+r_s*a**2*np.sin(theta)**2*(r**2-a**2*np.cos(theta)**2)) \
                             if not Sigma == 0 else np.inf

      Gamma_theta__t_t         = -C**2*r_s*a**2*np.sin(theta)*np.cos(theta)/Sigma**3 \
                                 if not Sigma == 0 else np.inf
      Gamma_theta__t_r         = 0
      Gamma_theta__t_theta     = 0
      Gamma_theta__t_phi       = C*r_s*a*r*(r**2+a**2)*np.sin(theta)*np.cos(theta)/Sigma**3 \
                                 if not Sigma == 0 else np.inf
      Gamma_theta__r_r         = a**2*np.sin(theta)*np.cos(theta)/(Sigma*Delta) \
                                 if not Sigma == 0 or Delta == 0 else np.inf
      Gamma_theta__r_theta     = r/Sigma if not Sigma == 0 else np.inf
      Gamma_theta__r_phi       = 0
      Gamma_theta__theta_theta = -a**2*np.sin(theta)*np.cos(theta)/Sigma if not Sigma == 0 else np.inf
      Gamma_theta__theta_phi   = 0
      Gamma_theta__phi_phi     = -np.sin(theta)*np.cos(theta)/Sigma**3*(A*Sigma+(r**2+a**2)*r_s*a**2*r*np.sin(theta)**2) \
                                 if not Sigma == 0 else np.inf

      Gamma_phi__t_t         = 0
      Gamma_phi__t_r         = C*r_s*a*(r**2-a**2*np.cos(theta)**2)/(2*Sigma**2*Delta) \
                               if not Sigma == 0 or Delta == 0 else np.inf
      Gamma_phi__t_theta     = -C*r_s*a*r*(np.cos(theta)/np.sin(theta))/Sigma**2 if not theta == 0 or Sigma==0 else np.inf
      Gamma_phi__t_phi       = 0
      Gamma_phi__r_r         = 0
      Gamma_phi__r_theta     = 0
      Gamma_phi__r_phi       = (2*r*Sigma**2+r_s*(a**4*np.sin(theta)**2*np.cos(theta)**2-r**2*(Sigma+r**2+a**2)))/(2*Sigma**2*Delta) \
                               if not Sigma == 0 or Delta == 0 else np.inf
      Gamma_phi__theta_theta = 0
      Gamma_phi__theta_phi   = (np.cos(theta)/np.sin(theta))/Sigma**2*(Sigma**2+r_s*a**2*r*np.sin(theta)**2) \
                               if not theta == 0 or Sigma == 0 else np.inf
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
