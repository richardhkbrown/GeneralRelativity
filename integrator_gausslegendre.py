import numpy as np

class Integrator_GaussLegendre:

   x_states = np.nan

   mF = []
   mX = []
   mNumStates = 0
   mDt = 0.01
   mThreshold = 1e-7
   mDamping = 1
   mMaxIterations = 100
   sq3 = np.sqrt(3)
   mA_Butcher = np.array([[1/4,1/4-sq3/6], \
                          [1/4+sq3/6,1/4]])
   mNumButcher = len( mA_Butcher[0,:] )
   mX_Butcher = []
   mK_Butcher = []
   mJ_Butcher = []
   mError = []
   mEta = 1e-8

   # define the state derivative function
   def DefineF( self, F ):
      self.mF = F
      return

   def SetX( self, X ):
      global mX,mNumStates
      if not len(X.shape):
         raise TypeError('X must be a vector')
      if not isinstance(X, np.ndarray):
         raise TypeError('X mist be a numpy vector')
      self.mX = X # same instance
      self.mNumStates = len( self.mX )
      return

   def SetDt( self, dT ):
      self.mDt = dT # same instance
      return
      
   def GuessK( self ):
      mF = self.mF
      mX = self.mX
      mNumButcher = self.mNumButcher
      mA_Butcher = self.mA_Butcher
      mDt = self.mDt
      mK_Butcher = self.mK_Butcher
      k_guess = mF( mX )
      x_guess = []
      for row in range( 0, mNumButcher ):
         x_guess.append( mX + np.sum(mA_Butcher[row,:])*mDt*k_guess )
      if len(mK_Butcher)==0:
         for row in range( 0, mNumButcher ):
            mK_Butcher.append( mF( x_guess[row] ) )
      else:
         for row in range( 0, mNumButcher ):
            mK_Butcher[row] =  mF( x_guess[row] )
      return

   def NextX( self ):
      mX = self.mX
      mDt = self.mDt
      mNumButcher = self.mNumButcher
      mA_Butcher = self.mA_Butcher
      mK_Butcher = self.mK_Butcher
      mX_Butcher = self.mX_Butcher
      for row in range( 0, mNumButcher ):
         newX = mX.copy()
         for col in range( 0, mNumButcher ):
            newX += mA_Butcher[row][col]*mK_Butcher[col]*mDt
         if len( mX_Butcher )<mNumButcher:
            mX_Butcher.append( newX )
         else:
            mX_Butcher[row] = newX
      return

   def Error( self ):
      mError = self.mError
      mNumButcher = self.mNumButcher
      mK_Butcher = self.mK_Butcher
      mF = self.mF
      mX_Butcher = self.mX_Butcher
      mError = []
      for row in range( 0, mNumButcher ):
         newK = mK_Butcher[row] - mF( mX_Butcher[row] )
         if len( mError ) == 0:
            mError = newK
         else:
            mError = np.hstack( (mError, newK) )
      mError = mError.reshape(-1,1)
      return

   def Jacobian( self, X, F, eta ):
      N = len( X )
      myJ = np.zeros([N,N])
      for col in range(0,N):
         myX = X.copy()
         myX[col] -= 0.5*eta
         myF0 = F(myX)
         myX = X.copy()
         myX[col] += 0.5*eta
         myF1 = F(myX)
         myF = (myF1 - myF0)/eta
         myJ[:,col] = myF.reshape(1,-1)
      return myJ

   def Jacobians( self ):
      global mJ_Butcher,mNumButcher,mX_Butcher,mF,mEta,mDt
      for row in range( 0, mNumButcher ):
         newJ = Jacobian( mX_Butcher[row], mF, mEta )
         if len( mJ_Butcher ) < mNumButcher:
            mJ_Butcher.append( newJ )
         else:
            mJ_Butcher[row] = newJ
      return

   def Step( self ):
      mError = self.mError
      mThreshold = self.mThreshold
      mJ_Butcher = self.mJ_Butcher
      mNumStates = self.mNumStates
      mDamping = self.mDamping
      mK_Butcher = self.mK_Butcher
      mMaxIterations = self.mMaxIterations
      mNumButcher = self.mNumButcher
      mDt = self.mDt
      mX = self.mX
      self.GuessK()
      self.NextX()
      self.Error()
      J_matrix = []
      J_row = []
      k_next = []
                  
      for iteration in range( 0, mMaxIterations ):
         if np.linalg.norm(mError) <= mThreshold:
            break
         self.NextX()
         self.Jacobians()
         for row in range( 0, mNumButcher ):
            for col in range( 0, mNumButcher ):
               if col==row:
                  J_row_new = np.eye( mNumStates ) \
                              -mDt*mA_Butcher[row][col]*mJ_Butcher[row]
               else:
                  J_row_new = -mDt*mA_Butcher[row][col]*mJ_Butcher[row]
               if len(J_row)==0 or  J_row.shape[1]<mNumButcher*mNumStates:
                  if len(J_row)==0:
                     J_row = J_row_new
                  else:
                     J_row = np.hstack( (J_row, J_row_new) )
               else:
                  J_row[0:mNumStates,col*mNumStates:(col+1)*mNumStates] = J_row_new
            if len(J_matrix)==0 or J_matrix.shape[0]<mNumButcher*mNumStates:
               if row==0:
                  J_matrix = J_row.copy()
               else:
                  J_matrix = np.vstack( (J_matrix,J_row) )
            else:
               J_matrix[row*mNumStates:(row+1)*mNumStates][:] = J_row
         for row in range( 0, mNumButcher ):
            if len( k_next )==0 or k_next.shape[0]<mNumButcher*mNumStates:
               if row==0:
                  k_next = mK_Butcher[row]
               else:
                  k_next = np.hstack( (k_next, mK_Butcher[row]) )
            else:
               k_next[row*mNumStates:(row+1)*mNumStates] = mK_Butcher[row]
         k_next = k_next.reshape(-1,1)
         try:
            k_next = k_next - \
                     mDamping * np.linalg.solve(J_matrix, mError)
         except Exception as e:
            print(e)
            breakpoint()
         k_next = k_next.reshape(1,-1)[0]
         for row in range( 0, mNumButcher ):
            mK_Butcher[row] = k_next[row*mNumStates:(row+1)*mNumStates]
         self.Error()
         
      if np.linalg.norm(mError) > mThreshold:
         raise TypeError('Newton did not converge by %d iterations.'%mMaxIterations)

      for row in range( 0, mNumButcher ):
         if row==0:
            sumK = mK_Butcher[row]
         else:
            sumK += mK_Butcher[row]
      mXpre = mX.copy()
      mX += (mDt / 2) * sumK
         
      return
