from numpy import *

class SolverExplicitRK(object):
    #---------------------------------------------------------------------------
    def __init__(self, cfl=None, deltaT=None, timeStepMode='', alfa={},
                 relTol=1e-4):

        if deltaT != None:
           self.__dt  = deltaT
           self.__cfl = 0.
        else:
           self.__cfl = cfl
           self.__dt = 0.

        if timeStepMode == 'local':
           self.__tModeLocal = 1
        elif timeStepMode == 'global':
           self.__tModeLocal = 0
        else:
           sys.exit(' *** ERROR: Unsupported time step mode. ***')
       
        self.__ntSub = len(alfa)
        self.__alfaSub = alfa
        self.__relTol = relTol

        self.__firstStep = 1
        self.__dw0 = None

        print " ... Solver parameters"
        print "     Explcit RK Solver:"
        print "        number of stages:", self.__ntSub
        print "        stage Coefficients:", self.__alfaSub
        print "        time step mode:         ", timeStepMode
        if deltaT != None:
            print "        time step size:         ", self.__dt
        else:
            print "        cfl:                    ", self.__cfl

    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    def getDeltaT(self,eqn,cfl):
        msh = eqn.domain()

        dt = cfl*msh.mass / eqn.spectralRadius()

        if not(self.__tModeLocal):
            #dt = min(dt)*ones(dt.shape)
            dt = min(dt)*ones(msh.nP)
        return dt
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    def solve(self,eqn,spaceSchm):
        msh = eqn.domain()

        cfl = self.__cfl
        if cfl != 0:
            #dt = self.__getDeltaT(cfl,tModeLocal,eqn.solution(),eqn,msh)
            dt = self.getDeltaT(eqn,cfl)
        else:
            dt = self.__dt* ones(msh.nP)
    
        eqn.storeSolution()
    
        # sub iterations
        for iSub in range(self.__ntSub):
    
            dw = eqn.computeResidual(spaceSchm)
            #for iV in range(eqn.nVar):
            #    print dw[iV,:]

            # update solution
            eqn.updateSolution(dw,(dt*self.__alfaSub[iSub]))
    
            # directly set values on boundaries
            eqn.correctBoundaryValues()
    
        #.......................................................................
        # end of subiterations
        #.......................................................................
    
        # return statistics
        wPrev = eqn.previousSolution()
        w = eqn.solution()
        if self.__firstStep:
           self.__dw0 = abs(wPrev-w).max(1)
           self.__firstStep = 0
        # compute max residual residuals and print to screen
        dwMax = abs((wPrev-w)).max(1)/self.__dw0
    
        if all(dwMax < self.__relTol):
            converged = 1
        else:
            converged = 0

        return converged, dwMax
