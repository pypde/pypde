from numpy import *

class SolverJamesonRK(object):
    convFluxList = ['ifsCentral2']
    diffFluxList = ['ifsCusp','ifsRoe']
    #---------------------------------------------------------------------------
    def __init__(self, cfl=None, deltaT=None, timeStepMode='', 
                 alfa={}, bta={},
                 relTol=1e-4, correctBCFaces=1, enforceStageBC=0):

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
       
        if len(alfa) != len(bta):
           sys.exit(' *** ERROR: Unsupported time step mode. ***')
 
        self.__ntSub = len(alfa)
        self.__alfaSub = alfa
        if bta[0] != 1.:
            print " *** Warning: First-stage diffusive coefficient set to 1"
            bta[0] = 1.
        self.__btaSub = bta
        self.__relTol = relTol
        self.__correctBCFaces = correctBCFaces
        self.__enforceStageBC = enforceStageBC

        self.__firstStep = 1
        self.__dw0 = None
        self.__fND = None

        print " ... Solver parameters"
        print "     Explcit RK Solver:"
        print "        number of stages", self.__ntSub
        print "        stage Coefficients:", self.__alfaSub
        print "        diffusive flux Coefficients:", self.__btaSub
        print "        time step mode:", timeStepMode
        if deltaT != None:
            print "        time step size:", self.__dt
        else:
            print "        cfl:", self.__cfl
        print "        relative tolerance:", self.__relTol
        print "        correct boundary face values:", self.__correctBCFaces
        print "        enforce BCs at each stage:", self.__enforceStageBC

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
        wPrev = eqn.previousSolution()
    

        # get interfaces solvers and split flux schemes into convective / diff
        ifsList = eqn.interfaceSolvers()
        convectiveFluxes = []
        diffusiveFluxes = []
        for ifs in ifsList:
            #if ifs.typeName in self.diffFluxList:
            if ifs.__class__.__name__ in self.diffFluxList:
                diffusiveFluxes.append(ifs)
            else:
                convectiveFluxes.append(ifs)
        
        if self.__firstStep:
            self.__fND = zeros([eqn.nVar,msh.nE])
        #else:
        #    self.__fND[:,:] = 0.

        # sub iterations
        for iSub in range(self.__ntSub):

            #dw = eqn.computeResidual(spaceSchm)

            w = eqn.solution()
            p = eqn.extraVariables()
            wL,wR,fL,fR = spaceSchm.spSol2FpFlx(eqn, w, p, msh)
            #wL,wR,fL,fR = spaceSchm.spSol2FpFlx(eqn, varDict, auxDict, msh)

            self.__fND *= (1. - self.__btaSub[iSub])
            fN = zeros([eqn.nVar,msh.nE])

            # interior faces
            for ifs in convectiveFluxes:
                fN += ifs.interfaceNormalFlux(wL,wR,fL,fR,
                                              msh.sF,msh,
                                              vstack([w,p]), \
                                              'interior')

            for ifs in diffusiveFluxes:
                self.__fND += self.__btaSub[iSub] \
                            * ifs.interfaceNormalFlux(wL,wR,fL,fR,
                                                      msh.sF,msh,
                                                      vstack([w,p]), \
                                                     'interior')

            fN += self.__fND
            dw = spaceSchm.fpFlx2SpRes(eqn,fN,msh)

            bcs = eqn.boundaryConditions()
            # boundary faces
            for iH in range(msh.numPatches()):
                dw[:,msh.facesOnPatch(bcs[iH].patchNum())] \
            += bcs[iH].faceNormalFlux()


            eqn.setSolution \
                (wPrev \
                 - dw*tile((dt*self.__alfaSub[iSub]) \
                           /msh.mass,[eqn.nVar,1]) \
                )

            # directly set values on boundaries
            
            if (self.__enforceStageBC and self.__correctBCFaces):
                eqn.correctBoundaryValues(itType='fractional')
    
        #.......................................................................
        # end of subiterations
        #.......................................................................
    
        if self.__correctBCFaces and (self.__enforceStageBC == 0):
            eqn.correctBoundaryValues(itType='full')

        # return statistics
        w = eqn.solution()
        if self.__firstStep:
           self.__dw0 = abs(wPrev-w).max(1)
           self.__firstStep = 0
        # compute max residual residuals and print to screen
        dwMax = abs((wPrev-w)).max(1)/self.__dw0
        eqn.storeMaxResiduals(dwMax)
    
        if all(dwMax < self.__relTol):
            converged = 1
        else:
            converged = 0

        return converged, dwMax
