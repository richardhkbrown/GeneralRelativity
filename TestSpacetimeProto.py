import numpy as np
from spacetime_minkowski import SpaceTime_Minkowski
from spacetime_schwarschild import SpaceTime_Schwarschild
from spacetime_kerr import SpaceTime_Kerr
import matplotlib.pyplot as plt
from integrator_gausslegendre import Integrator_GaussLegendre
import sys

metric = "Schwarschild"
maxView = 50
h0 = 0.1
               
r2o2 = np.sqrt(2)/2

maxTime = 249
inits = []
inits.append({"x":30,"y":0,"z":0,"hasMass":True, \
              "v_x":0,"v_y":0.578,"v_z":0,"h":0.1})
inits.append({"x":0,"y":-30,"z":0,"hasMass":True, \
              "v_x":0.578,"v_y":0,"v_z":0,"h":h0})
inits.append({"x":-30,"y":0,"z":0,"hasMass":True, \
              "v_x":0,"v_y":-0.578,"v_z":0,"h":h0})
inits.append({"x":0,"y":30,"z":0,"hasMass":True, \
              "v_x":-0.578,"v_y":0,"v_z":0,"h":h0})

fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True)
ax1.set_aspect('equal')
ax1.grid(True)
ax1.set_xlim(-3*maxView-1,3*maxView+1)
ax1.set_ylim(-3*maxView-1,3*maxView+1)
ax2.set_aspect('equal')
ax2.grid(True)

# plot light
xdata = np.array([])
ydata = np.array([])
xdataLorentz = np.array([])
ydataLorentz = np.array([])
line1, = ax1.plot( [], [], linewidth=0.5 )
line2, = ax2.plot( [], [], linestyle='None', markersize=10 )

figSize = fig.get_size_inches()
fig.set_size_inches( [1.5*figSize[0],figSize[1]] )

plt.show(block=False)
plt.ion()
 
handles = []
st_temp = eval("SpaceTime_"+metric+"()")
for iobject, init in enumerate( inits ):
   st_temp.x_states = np.nan # re-init temporary metric
   x = init["x"]
   y = init["y"]
   z = init["z"]
   v_x = init["v_x"]
   v_y = init["v_y"]
   v_z = init["v_z"]
   hasMass = init["hasMass"]
   h = init["h"]
   v_mag = np.sqrt( v_x**2 + v_y**2 + v_z**2 )
   x_temp = st_temp.CartesianToStates( 0, v_x, v_y, v_z, \
                                       0, x, y, z )
   v_t, v_r, v_theta, v_phi, \
      t, r, theta, phi = st_temp.StatesToVariables( x_temp )
   delta = 1
   for icount in range(1000):
      v_t = np.nan
      X = st_temp.VariablesToStates( v_t, v_r, v_theta, v_phi, \
                                     t, r, theta, phi, hasMass )
      if not hasMass:
         v_achieved = v_mag
         break
      v_t,v_x,v_y,v_z,t,x,y,z = \
         st_temp.StatesToCartesian( X )
      v_achieved = np.sqrt( v_x**2 + v_y**2 + v_z**2 ) / v_t    
      if v_achieved < v_mag:
         scl = 1 + delta
      else:
         scl = 1 - delta
      v_r *= scl
      v_theta *= scl
      v_phi *= scl
      delta *= 0.9
   print("v_achieved",v_achieved)
   lineInertial, = ax1.plot( [], [], linewidth=0.5 )
   line, = ax2.plot( [], [], linestyle='', marker='o', markersize=10 )
   handles.append( {"metric":eval("SpaceTime_"+metric+"()"), \
                    "integrator":Integrator_GaussLegendre(), \
                    "states":np.nan, \
                    "worldlineGodsEye":np.array([]), \
                    "lineGodsEye":lineInertial, \
                    "lineLocal":line} )
   handles[-1]["metric"].hasMass = hasMass
   handles[-1]["metric"].x_states = X.copy()
   handles[-1]["metric"].name = "object%2d"%iobject
   handles[-1]["integrator"].DefineF( handles[-1]["metric"].F )
   handles[-1]["integrator"].SetX( handles[-1]["metric"].x_states )
   handles[-1]["integrator"].SetDt( h )
   vt,vx,vy,vz,t,x,y,z = handles[-1]["metric"].StatesToCartesian()
   worldline = np.array( [[vt,vx,vy,vz,t,x,y,z]] )
   worldlineLorentz = handles[0]['metric'].LorentzTransform( worldline )
   handles[-1]["worldlineGodsEye"] = worldline

# initialize all instances
globalMetric = eval("SpaceTime_"+metric)
globalMetric.C = 1
if "Schwarschild" == metric or "Kerr" == metric:
   G = 1
   M = 10
   globalMetric.G = G
   globalMetric.M = M
   globalMetric.r_s = 2*G*M/st_temp.C**2
if metric == "Kerr":
   J = 10
   globalMetric.J = J
   globalMetric.a = J/(st_temp.M*st_temp.C)

# plot schwartzchild radius
if 'r_s' in dir(handles[0]["metric"]):
   r_s = handles[0]["metric"].r_s
else:
   r_s = 1.
x_r_s = []
y_r_s = []
for ang in range(0,361):
   x_r_s.append( r_s*np.cos(np.deg2rad(ang)))
   y_r_s.append( r_s*np.sin(np.deg2rad(ang)))
ax1.plot( x_r_s, y_r_s, linewidth=0.5 )

# Integrate
dt = 0.1
for iLoop in range(100000):
   time = dt*iLoop
   if time > maxTime:
      break
   for handle in handles:

      if r_s==r_s:
         if handle["metric"].x_states[5] <= r_s:
            continue
      
      while handle["metric"].x_states[4] < time:
         handle["integrator"].Step()
         handle["metric"].NormalizeVelocity()
      vt,vx,vy,vz,t,x,y,z = handle["metric"].StatesToCartesian()
      worldline = np.array( [[vt,vx,vy,vz,t,x,y,z]] )
      handle['worldlineGodsEye'] = np.vstack( (handle['worldlineGodsEye'], \
                                               worldline[0]) )
      if iLoop%100==0:
         handle['lineGodsEye'].set_data( handle['worldlineGodsEye'][:,5], \
                                         handle['worldlineGodsEye'][:,6] )
         plt.pause(0.01)

# Plot
for handle in handles:
   st_temp.PlotWorldline( handle['lineGodsEye'], handle['worldlineGodsEye'] )
         
# setup light integrator
st_temp.hasMass = False
int_temp = Integrator_GaussLegendre()
int_temp.DefineF( st_temp.F )
int_temp.SetX( st_temp.x_states )
int_temp.SetDt( h )

def Observation( handle0, handle ):

   global maxView, int_temp, line1

   # create light world lines around prime world line
   primeWorldLine = handle0['worldlineGodsEye'][-1]
   v_t0 = primeWorldLine[0]
   v_x0 = primeWorldLine[1]
   v_y0 = primeWorldLine[2]
   v_z0 = primeWorldLine[3]
   t0 = primeWorldLine[4]
   x0 = primeWorldLine[5]
   y0 = primeWorldLine[6]
   z0 = primeWorldLine[7]
   lightWorldLines = {}
   lightWorldLinesLorentz = {}
   r_s = 1.
   if 'r_s' in dir(st_temp):
      r_s = st_temp.r_s

   def GetDistanceAll( angle, handle0, handle1 ):

      global maxView, int_temp, line1, xdata, ydata, \
             xdataLorentz, ydataLorentz
      
      primeWorldLine = handle0['worldlineGodsEye'][-1]
      v_t0 = primeWorldLine[0]
      v_x0 = primeWorldLine[1]
      v_y0 = primeWorldLine[2]
      v_z0 = primeWorldLine[3]
      t0 = primeWorldLine[4]
      x0 = primeWorldLine[5]
      y0 = primeWorldLine[6]
      z0 = primeWorldLine[7]
      try:
         x_temp = st_temp.CartesianToStates( \
                  1,np.cos(angle),np.sin(angle),0,\
                  t0,x0,y0,z0 )
      except Exception as e:
         print(e)
         breakpoint()
      v_t, v_r, v_theta, v_phi, \
         t, r, theta, phi = st_temp.StatesToVariables( x_temp )
      v_t = np.nan
      st_temp.VariablesToStates( v_t, v_r, v_theta, v_phi, \
                                 t, r, theta, phi )
      v_t,v_x,v_y,v_z,t,x,y,z = st_temp.StatesToCartesian()
      #lightWorldLine = np.array([[v_t,v_x,v_y,v_z,t,x,y,z]])
      lightWorldLine = []
      for iLoop in range(100000):
         mK_Butcher = int_temp.mK_Butcher
         mNumButcher = int_temp.mNumButcher
         if len(mK_Butcher)>0:
            for row in range( 0, mNumButcher ):
               if not all(mK_Butcher[row]==mK_Butcher[row]):
                  print("A mK_Butcher[row]",mK_Butcher[row])
                  breakpoint()
         if any(np.isinf(X)):
            print("Ab",int_temp.mX)
            breakpoint()
         int_temp.Step()
         if (np.abs(int_temp.mX[5])<=r_s):
            break
         st_temp.NormalizeVelocity()
         if int_temp.mX[4]>t0+10.:
            v_t,v_x,v_y,v_z,t,x,y,z = st_temp.StatesToCartesian()
            newValue = np.array( [[vt,vx,vy,vz,t,x,y,z]] )
            if len(lightWorldLine)==0:
               lightWorldLine = newValue
            else:
               lightWorldLine = np.vstack( (lightWorldLine, newValue[0]) )
         if x>maxView or x<-maxView or \
            y>maxView or y<-maxView:
            break

      objWorldLine = handle1['worldlineGodsEye']
      lightpos = lightWorldLine[:,4:8]
      lightpos[:,0] = 2.*lightpos[0,0] - lightpos[:,0]
      objpos = objWorldLine[:,4:8]
      deltas = objpos[:,None] - lightpos
      s_sq = (st_temp.C*deltas[:,:,0])**2 + \
             (deltas[:,:,1]**2+deltas[:,:,2]**2+deltas[:,:,3]**2)
      differences = np.abs(s_sq)
      min_index_flat = np.argmin(differences)
      row, col = np.unravel_index(min_index_flat, differences.shape)
      objPoint = objWorldLine[row,:]
      lightPoint = lightWorldLine[col,:]

      st_temp.PlotWorldline( line1, lightWorldLine )      
      plt.pause(0.01)

      return np.min(differences), objPoint, lightPoint

   def GetDistance( angle, handle0, handle1 ):
      distance,_,_ = GetDistanceAll( angle, handle0, handle1 )
      return distance

   # CUSTOM!!!

   if False:
      angle = 235*np.pi/180.
      thisDistance = GetDistance(angle,handle0,handle)
      plt.ion()
      plt.pause(0.5)
      print("thisDistance %9.3f"%thisDistance)
      sys.exit()

   # CUSTOM!!!

   x = np.array([])
   f = np.array([])
   for angle0 in np.arange(0., 361.*np.pi/180., 1.*np.pi/180. ):
      angle = angle0
      x = np.append(x,angle)
      f = np.append(f,GetDistance(angle,handle0,handle))
   sidx = np.argsort(f)
   x = x[sidx]
   f = f[sidx]
   x0 = x[0]
   f0 = GetDistance(x0,handle0,handle)

   x = np.array([])
   f = np.array([])
   for angle0 in np.arange(0., 2.0*np.pi/180., 0.1*np.pi/180. ):
      angle = x0 + angle0
      x = np.append(x,angle)
      f = np.append(f,GetDistance(angle,handle0,handle))      
      angle = x0 - angle0
      x = np.append(x,angle)
      f = np.append(f,GetDistance(angle,handle0,handle))
   sidx = np.argsort(f)
   x = x[sidx]
   f = f[sidx]
   x0 = x[0]
   f0 = GetDistance(x0,handle0,handle)

   x = np.array([])
   f = np.array([])
   for angle0 in np.arange(0., 0.2*np.pi/180., 0.01*np.pi/180. ):
      angle = x0 + angle0
      x = np.append(x,angle)
      f = np.append(f,GetDistance(angle,handle0,handle))      
      angle = x0 - angle0
      x = np.append(x,angle)
      f = np.append(f,GetDistance(angle,handle0,handle))
   sidx = np.argsort(f)
   x = x[sidx]
   f = f[sidx]
   x0 = x[0]
   f0 = GetDistance(x0,handle0,handle)

   print("%9.3f"%f0)

   minAngle = x[0]
   return GetDistanceAll(minAngle,handle0,handle)

for handle in handles:
   distance, objPoint, lightPoint = Observation( handles[0], handle )
   lightPoint = np.array([lightPoint])
   lightWorldLineLorentz = handles[0]['metric'].LorentzTransform( lightPoint )
   handle['lineLocal'].set_data( [lightWorldLineLorentz[0][5]],
                                 [lightWorldLineLorentz[0][6]] )
   print("%9.3f %9.3f"%(lightWorldLineLorentz[0][5],lightWorldLineLorentz[0][6]))
   plt.ion()
   plt.pause(0.1)
