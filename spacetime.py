import numpy as np

class SpaceTime:

   x_states = np.nan
   hasMass = True
   name = "unknown"

   # v_t     = (t_1 - t_0) / lambda
   # v_r     = (r_1 - r_0) / lambda
   # v_theta = (theta_1 - theta_0) / lambda
   # v_phi   = (phi_1 - phi_0) / lambda
   # t       = t_1
   # r       = r_1
   # theta   = theta_1
   # phi     = phi_1
   #
   # for mass-like objects, lamda is tau (proper time) or the time experienced by the object
   # vs t (coordinate time) as observed from in observer an infinite distance away in flat-space

   def StatesToVariables( self, X = np.nan ):
      if np.isnan(X).any():
         X = self.x_states
      v_t     = X[0]
      v_r     = X[1]
      v_theta = X[2]
      v_phi   = X[3]
      t       = X[4]
      r       = X[5]
      theta   = X[6]
      phi     = X[7]
      return v_t, v_r, v_theta, v_phi, \
             t, r, theta, phi

   def VariablesToStates( self, v_t, v_r, v_theta, v_phi, \
                          t, r, theta, phi, hasMass = np.nan ):
      if hasMass==hasMass:
         self.hasMass = hasMass
      x_states = self.x_states
      if not v_t==v_t:
         C = self.C
         # these simplifications can be made only if coordinates are orthogonal
         # ds^2 = -c^2*dtau^2 = g__dt_dt*dt^2 + g__dr_dr*dr^2 + g__dphi_dphi*dphi^2 + g__dtheta_dtheta*dtheta^2
         # -c^2 = g__dt_dt*dt^2/dtau^2 + g__dr_dr*dr^2/dtau^2 + g__dphi_dphi*dphi^2/dtau^2 + g__dtheta_dtheta*dtheta^2/dtau^2
         # v_t = dt/dtau
         # v_r = dr/dtau
         # v_phi = dphi/dtau
         # v_theta = dtheta/dtau
         # -c^2 = g__dt_dt*v_t^2 + g__dr_dr*v_r^2 + g__dphi_dphi*v_phi^2 + g__dtheta_dtheta*v_theta^2
         g__dt_dt,g__dt_dr,g__dt_dtheta,g__dt_dphi,g__dr_dr, \
            g__dt_dtheta,g__dr_dphi,g__dtheta_dtheta,g__dtheta_dphi, \
            g__dphi_dphi = self.MetricTensor( r, theta, phi )
         if self.hasMass:
            # v_t^2 = -c^2 - g__dr_dr*v_r^2 - g__dphi_dphi*v_phi^2 - g__dtheta_dtheta*v_theta^2
            v_t_2 = (-C**2 - g__dr_dr*v_r**2 - g__dphi_dphi*v_phi**2 \
                                 - g__dtheta_dtheta*v_theta**2) / g__dt_dt
            v_t = np.sqrt( v_t_2 )
         else:
            v_t = 0
            # ds^2 = 0
            # 0 = g__dt_dt*dt^2/dlamda^2 + g__dr_dr*dr^2/dlamda^2 + g__dphi_dphi*dphi^2/dlamda^2 + g__dtheta_dtheta*dtheta^2/dlamda^2
            v_t_2 = ( -g__dr_dr*v_r**2 - g__dphi_dphi*v_phi**2 \
                      - g__dtheta_dtheta*v_theta**2 ) / g__dt_dt
            if not v_t_2==0:
               # assert d_lamda = d_t, so v_t = 1
               v_t = np.sqrt( v_t_2 )
               v_r /= v_t
               v_theta /= v_t
               v_phi /= v_t
               v_t /= v_t
            else:
               print("v_t_2",v_t_2,"v_r",v_r,"v_phi",v_phi,"v_theta",v_theta)
               raise ValueError("Time velocity is zero")
      if np.isnan(x_states).any():
         x_states = np.array([v_t, v_r, v_theta, v_phi, \
                              t, r, theta, phi])
         x_states = x_states.reshape(1,-1)[0]
      else:
         x_states[:] = [v_t, v_r, v_theta, v_phi, \
                        t, r, theta, phi]

      if np.isscalar( self.x_states ):
         self.x_states = x_states
      return x_states

   def F( self, X ):
      v_t, v_r, v_theta, v_phi, \
             t, r, theta, phi = self.StatesToVariables( X )
      v_v = np.array([v_t*v_t      ,v_r*v_t     ,v_theta*v_t     ,v_phi*v_t    , \
                      v_t*v_r      ,v_r*v_r     ,v_theta*v_r     ,v_phi*v_r    , \
                      v_t*v_theta  ,v_r*v_theta ,v_theta*v_theta ,v_phi*v_theta, \
                      v_t*v_phi    ,v_r*v_phi   ,v_theta*v_phi   ,v_phi*v_phi  ])

      Gamma_t,Gamma_r,Gamma_theta,Gamma_phi = self.ChristoffelSymbols(r,theta,phi)
      # The derivative of the states + ChristoffelSymbols*partial derivatives of states = 0
      # The derivative of the states = -ChristoffelSymbols*partial derivatives of states
      # The derivatives are with respect to dLambda
      dvt_dlamda = np.nan_to_num( -Gamma_t ) @ v_v
      dvr_dlamda = np.nan_to_num( -Gamma_r ) @ v_v
      dvtheta_dlamda = np.nan_to_num( -Gamma_theta ) @ v_v
      dvphi_dlamda = np.nan_to_num( -Gamma_phi ) @ v_v
      if np.isnan(dvt_dlamda):
         dvt_dlamda = 0.0
         for idx, a in enumerate(np.nan_to_num( -Gamma_t )):
            if a==0.0 and np.isinf(v_v[idx]):
               dvt_dlamda += 0.0
            else:
               dvt_dlamda += a*v_v[idx]
      if np.isnan(dvr_dlamda):
         dvr_dlamda = 0.0
         for idx, a in enumerate(np.nan_to_num( -Gamma_r )):
            if a==0.0 and np.isinf(v_v[idx]):
               dvr_dlamda += 0.0
            else:
               dvr_dlamda += a*v_v[idx]
      if np.isnan(dvtheta_dlamda):
         dvtheta_dlamda = 0.0
         for idx, a in enumerate(np.nan_to_num( -Gamma_theta )):
            if a==0.0 and np.isinf(v_v[idx]):
               dvtheta_dlamda += 0.0
            else:
               dvtheta_dlamda += a*v_v[idx]
      if np.isnan(dvphi_dlamda):
         dvphi_dlamda = 0.0
         for idx, a in enumerate(np.nan_to_num( -Gamma_phi )):
            if a==0.0 and np.isinf(v_v[idx]):
               dvphi_dlamda += 0.0
            else:
               dvphi_dlamda += a*v_v[idx]
               
      k = np.array([dvt_dlamda,dvr_dlamda,dvtheta_dlamda,dvphi_dlamda,v_t,v_r,v_theta,v_phi])

      return k

   def StatesToCartesian( self, X = np.nan ):
      v_t, v_r, v_theta, v_phi, \
         t, r, theta, phi = self.StatesToVariables( X )
      r_dot = v_r
      theta_dot = v_theta
      phi_dot = v_phi
      x = r*np.sin(theta)*np.cos(phi)
      y = r*np.sin(theta)*np.sin(phi)
      z = r*np.cos(theta)
      v_x = r_dot*np.sin(theta)*np.cos(phi)+r*theta_dot*np.cos(theta)*np.cos(phi)-r*phi_dot*np.sin(theta)*np.sin(phi)
      v_y = r_dot*np.sin(theta)*np.sin(phi)+r*theta_dot*np.cos(theta)*np.sin(phi)+r*phi_dot*np.sin(theta)*np.cos(phi)
      v_z = r_dot*np.cos(theta)-r*theta_dot*np.sin(theta)
      return v_t,v_x,v_y,v_z,t,x,y,z

   def CartesianToStates( self, v_t, v_x, v_y, v_z, \
                          t, x, y, z ):
      r = np.sqrt( x**2 + y**2 + z**2 )
      if not r == 0:
         theta = np.acos( z / r )
      else:
         theta = 0.
      phi = np.atan2( y, x )
      if not r == 0:
         v_r = ( x*v_x + y*v_y + z*v_z ) / r
         v_theta = ( x*z*v_x + y*z*v_y - (x**2 + y**2)*v_z ) / \
                   ( r**2 * np.sqrt(x**2 + y**2) )
         v_phi = ( x*v_y - y*v_x ) / \
                 ( r * np.sin( theta) * np.sqrt( x**2 + y**2 ) )         
      elif not ( v_x**2 + v_y**2 + v_z**2 ) == 0:
         v_r = np.sqrt( v_x**2 + v_y**2 + v_z**2 )
         v_theta = 0
         v_phi = 0
         phi = np.atan2( v_y, v_x )
         theta = np.acos( v_z / v_r )
      else:
         v_r = 0
         v_theta = 0
         v_phi = 0

      return [v_t, v_r, v_theta, v_phi, t, r, theta, phi]

   def NormalizeVelocity( self, X = np.nan ):
      if np.isnan(X).any():
         X = self.x_states         
      v_t, v_r, v_theta, v_phi, \
         t, r, theta, phi = self.StatesToVariables( X )
      g__dt_dt,g__dt_dr,g__dt_dtheta,g__dt_dphi,g__dr_dr, \
            g__dt_dtheta,g__dr_dphi,g__dtheta_dtheta,g__dtheta_dphi, \
            g__dphi_dphi = self.MetricTensor( r, theta, phi )      
      nc2 = g__dt_dt*v_t**2+g__dr_dr*v_r**2+g__dphi_dphi*v_phi**2+g__dtheta_dtheta*v_theta**2
      if self.hasMass:
         # -c**2 = g__dt_dt*v_t**2+g__dr_dr*v_r**2+g__dphi_dphi*v_phi**2+g__dtheta_dtheta*v_theta**2
         # -c**2 - g__dt_dt*v_t**2 = K*(g__dr_dr*v_r**2+g__dphi_dphi*v_phi**2+g__dtheta_dtheta*v_theta**2)
         # K = -(c**2 + g__dt_dt*v_t**2) / (g__dr_dr*v_r**2+g__dphi_dphi*v_phi**2+g__dtheta_dtheta*v_theta**2)
         den = g__dr_dr*v_r**2+g__dphi_dphi*v_phi**2+g__dtheta_dtheta*v_theta**2
         if not den==0:
            K = -(self.C**2 + g__dt_dt*v_t**2) / den
         else:
            K = 1
      else:
         # 0 = g__dt_dt*v_t**2+g__dr_dr*v_r**2+g__dphi_dphi*v_phi**2+g__dtheta_dtheta*v_theta**2
         # -g__dt_dt*v_t**2 = K*(g__dr_dr*v_r**2+g__dphi_dphi*v_phi**2+g__dtheta_dtheta*v_theta**2)
         # K = -g__dt_dt*v_t**2 / (g__dr_dr*v_r**2+g__dphi_dphi*v_phi**2+g__dtheta_dtheta*v_theta**2)
         den = g__dr_dr*v_r**2+g__dphi_dphi*v_phi**2+g__dtheta_dtheta*v_theta**2
         if not den==0:
            K = -g__dt_dt*v_t**2 / den
         else:
            K = 1
      if K>=0:
         K1 = np.sqrt(K)
         v_r *= K1
         v_theta *= K1
         v_phi *= K1         
      X[1] = v_r
      X[2] = v_theta
      X[3] = v_phi

   def LorentzTransform( self, otherWorldline ):
      #v_t0,v_x0,v_y0,v_z0,t0,x0,y0,z0 = self.StatesToCartesian( self.x_states )
      v_t0,v_x0,v_y0,v_z0,t0,x0,y0,z0 = self.StatesToCartesian()
      v = np.array([v_x0,v_y0,v_z0]) / v_t0
      V = np.linalg.norm( v )
      gamma = 1 / np.sqrt( 1 - V**2/self.C**2 ) \
              if not V**2/self.C**2 >= 1 else 1.0
      r = otherWorldline[:,5:8] - np.array([x0,y0,z0])
      t = otherWorldline[:,4] - t0
      t_prime = gamma*(t-(r@v)/self.C**2)
      tempV = ((gamma-1)*(r@v)/V**2-gamma*t)
      tempV = tempV.reshape(-1, 1)
      r_prime = r+tempV*v
      u = otherWorldline[:,1:4] / otherWorldline[:,0].reshape(-1,1)
      u_prime = 1 / ( 1 - ( u @ v ).reshape(-1,1) / self.C**2 ) * \
                ( u / gamma - v + 1 / self.C**2 * \
                gamma / ( gamma + 1 ) * \
                ( u @ v ).reshape(-1,1) * v )
      U_PRIME = np.linalg.norm( u_prime , axis=1, keepdims=True )
      vrat = U_PRIME**2 /  self.C**2
      vrat[ vrat>1 ] = 1
      vt_prime = np.sqrt( 1 - vrat )
      return np.hstack((np.hstack((vt_prime.reshape(-1,1),u_prime)),
                        np.hstack((t_prime.reshape(-1,1),r_prime))))

   def PlotWorldline( self, line, worldline ):
      line.set_data( worldline[:,5], worldline[:,6] )
      ax = line.axes
      for dt in [0.1,1,2,5,10,25,50,100,250,500,1000]:
         textTimesPre = dt*np.round(worldline[:,4]/dt)
         idxs = np.unique(textTimesPre, return_index=True)[1]
         textTimes = [textTimesPre[idx] for idx in sorted(idxs)]
         if len(textTimes)>1:
            if worldline[-1,4] > worldline[0,4]:
               x_a = np.interp( textTimes[0], worldline[:,4], worldline[:,5] )
               y_a = np.interp( textTimes[0], worldline[:,4], worldline[:,6] )
               x_b = np.interp( textTimes[1], worldline[:,4], worldline[:,5] )
               y_b = np.interp( textTimes[1], worldline[:,4], worldline[:,6] )
            else:
               x_a = np.interp( textTimes[0], worldline[::-1,4], worldline[::-1,5] )
               y_a = np.interp( textTimes[0], worldline[::-1,4], worldline[::-1,6] )
               x_b = np.interp( textTimes[1], worldline[::-1,4], worldline[::-1,5] )
               y_b = np.interp( textTimes[1], worldline[::-1,4], worldline[::-1,6] )
            deltaP = np.sqrt((x_a-x_b)**2+(y_a-y_b)**2)
            if deltaP > 2:
               break
      color = line.get_color()
      if 'annotations' in dir(line):
         for annotation in line.annotations :
            annotation.remove()
      line.annotations = []  
      for textTime in textTimes:
         if worldline[-1,4] > worldline[0,4]:
            xt = np.interp( textTime, worldline[:,4], worldline[:,5] )
            yt = np.interp( textTime, worldline[:,4], worldline[:,6] )
         else:
            xt = np.interp( textTime, worldline[::-1,4], worldline[::-1,5] )
            yt = np.interp( textTime, worldline[::-1,4], worldline[::-1,6] )
         annotation = ax.annotate( "%.0f"%textTime, xy=(xt,yt), xytext=(0, 0),
                      textcoords='offset points', fontsize=6,
                      color=color, ha='center', va='center' )
         line.annotations.append( annotation )
