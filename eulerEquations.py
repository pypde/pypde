from bcConstructor import *
from refStateConstructor import *
from ifsConstructor import *
from numpy import *


class EulerEquations(object):
    nD = 2
    nVar = 4
    type = 'Euler'
    #---------------------------------------------------------------------------
    # constructor
    #
    def __init__(self,msh,modelDict,refStateDict,bcDict,ifsDict):
        #self.__gasModel = constructThermoModel(modelDict)
        #self.__gma = self.__gasModel.gma(w?)
        self.__gma = modelDict.get('gamma')
        self.__refState = constructRefState(self,refStateDict)
        #self.__interfaceSolvers, self.__interfaceSolversArgs \
        self.__interfaceSolvers = constructInterfaceSolver(self,msh,ifsDict)
        self.__boundaryConditions = constructBC(self,self.__refState,msh,bcDict)
        self.__w = empty([self.nVar,msh.nP])
        self.__p = empty([1,msh.nP])
        self.__wPrev = None
        self.__msh = msh
        self.__maxRes = .5/(finfo(double).eps)*ones([self.nVar,1])
 
        self.varDict = {}
        self.varDict['rho']= self.rho
        self.varDict['rhou'] = self.rhou
        self.varDict['rhov'] = self.rhov
        self.varDict['rhoe'] = self.rhoe
        self.varDict['p'] = self.p
        self.varDict['mach'] = self.mach

        self.varOrder = {}
        self.varOrder['rho'] = 0
        self.varOrder['rhou'] = 1
        self.varOrder['rhov'] = 2
        self.varOrder['rhoe'] = 3
        self.varAuxOrder = {}
        self.varAuxOrder['p'] = 0
        self.varAuxOrder['mach'] = 1
        self.varAuxOrder['rhoh'] = 2

        print " ... Model parameters"
        print "     Equation:"
        print "         Euler equations"
        print "     Gas:"
        print "         ratio of specific heats:", self.__gma

    #---------------------------------------------------------------------------
    # access functions
    #---------------------------------------------------------------------------
    def get(self,fieldName):
        return self.varDict.get(fieldName)
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    def domain(self):
        return self.__msh
    #---------------------------------------------------------------------------
    def gma(self):
        return self.__gma
    #---------------------------------------------------------------------------
    def refState(self):
        return self.__refState
    #---------------------------------------------------------------------------
    def rho(self):
        return self.__w[0,:]
    #---------------------------------------------------------------------------
    def rhou(self):
        return self.__w[1,:]
    #---------------------------------------------------------------------------
    def rhov(self):
        return self.__w[2,:]
    #---------------------------------------------------------------------------
    def rhoe(self):
        return self.__w[3,:]
    #---------------------------------------------------------------------------
    def p(self):
        return self.__p[0,:]
    #---------------------------------------------------------------------------
    def h(self):
        return (self.__w[3,:] + self.__p[0,:])/(self.__w[0,:])
    #---------------------------------------------------------------------------
    def solution(self):
        return self.__w
    #---------------------------------------------------------------------------
    def maxResiduals(self):
        return self.__maxRes
    #---------------------------------------------------------------------------
    def previousSolution(self):
        return self.__wPrev
    #---------------------------------------------------------------------------
    def extraVariables(self):
        return self.__p
    #---------------------------------------------------------------------------
    def allVariables(self):
        return vstack([self.__w,self.__p])
    #---------------------------------------------------------------------------
    def boundaryConditions(self):
        return self.__boundaryConditions
    #---------------------------------------------------------------------------
    def interfaceSolvers(self):
        return self.__interfaceSolvers
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    # calculation functions
    #---------------------------------------------------------------------------
    def updateExtraVariables(self, w=None):
        p = self.updateP(w)
        return p
    #---------------------------------------------------------------------------
    #def updateP(self,w):
    def updateP(self,w=None):
        if w == None:
            w = self.__w
            p = self.__p
        else:
            p = zeros([1,w.shape[1]])

        p[0,:] = (self.__gma-1) * (w[3,:] - .5*(w[1,:]**2 + w[2,:]**2)/w[0,:])
        if any(p[0,:]<0):
            print ' *** WARNING: negative pressure set to zero in calculation ***'
            p[p<0] = 0.

        return p[0,:]
    #---------------------------------------------------------------------------
    #def mach(self, w):
    def mach(self):
        w = self.__w
        return sqrt((w[1,:]**2 + w[2,:]**2) / (self.__gma*w[0,:]*self.__p[0,:]))
    #---------------------------------------------------------------------------
    def storeSolution(self):
        self.__wPrev = self.__w.copy()
    #---------------------------------------------------------------------------
    def storeMaxResiduals(self,res):
        self.__maxRes = res.copy()
    #---------------------------------------------------------------------------
    #def initSolution(self, msh, refState=None):
    def initSolution(self, refState=None):
        if refState == None:
            refState = self.__refState
        rho = refState.rho()
        p = refState.p()
        w = self.__w
        w[0,:] = refState.rho()
        w[1,:] = refState.u()
        w[2,:] = refState.v()
        w[3,:] = rho * refState.h() - p
        #w[4,:] = p
        self.__p[0,:] = p
        print '     density:               ', w[0,0]
        print '     x-momentum:            ', w[1,0]
        print '     y-momentum:            ', w[2,0]
        print '     energy per unit volume:', w[3,0]
        print '     pressure:              ', self.__p[0,0]
        print ' done.'

        # return w
    #---------------------------------------------------------------------------
    def setSolution(self, wNew):

        self.__w = wNew.copy()

        if any(self.__w[0,:] < 0.):
            print ' *** WARNING: negative density set to zero during update ***'
            self.__w[0,(self.__w[0,:]<0.)] = 0.

        self.updateP()
    #---------------------------------------------------------------------------
    def pointFlux(self,w,p):
    #def pointFlux(self):

        #w = self.__w
        #p = self.__p
        #nP = self.__msh.nP
        nP = w.shape[1]

        #p = w[4,:]
        #p = wAux
        rho = w[0,:]
        rhou = w[1,:]
        rhov = w[2,:]
        rhoePp = w[3,:] + p

        u = w[1,:]/rho
        v = w[2,:]/rho

        f = zeros([self.nD,self.nVar,nP])
        f[0,0,:] = rhou
        f[0,1,:] = rhou*u + p
        f[0,2,:] = rhou*v
        f[0,3,:] = (rhoePp)*u

        f[1,0,:] = rhov
        f[1,1,:] = rhou*v
        f[1,2,:] = rhov*v + p
        f[1,3,:] = (rhoePp)*v

        return f
    #---------------------------------------------------------------------------
    #def spectralRadius(self,w,msh):
    def spectralRadius(self):
        # sF: nDim x nF array of face areas
        # average face-normal component of (u+c) using central difference

        msh = self.__msh
        e0 = msh.e0
        e1 = msh.e1
        sF = msh.sF

        w = self.__w
        p = self.__p

        uDotN  = sum( sF \
                     *( (w[1:3,e0] / tile(w[0,e0],[2,1])) \
                       +(w[1:3,e1] / tile(w[0,e1],[2,1]))),
                     0)

        cF = sqrt(self.__gma*sum(sF**2,0)*(  p[0,e0]/w[0,e0] \
                                               + p[0,e1]/w[0,e1]))
        radF = abs(uDotN) + cF
        rad = array([sum(radF[hstack([msh.egClouds0[ii],msh.egClouds1[ii]])]) \
                     for ii in range(msh.nP)])

        return rad
    #---------------------------------------------------------------------------
    #def computeDerivativeMatrix(self,w,p,spaceSchm,msh):

    #    move to SOLVER layer (Jac-Free-Implicit!!)
    #    epsMSqrt = sqrt(finfo(double).eps)
    #    delta = 1e-2
 
    #    store = True
    #    wL,wR,fL,fR = spaceSchm.spSol2FpFlx(self, w, p, msh, store) # repeated...
    #    #dwL,dwR,dfL,dfR = spaceSchm.getGrad(self, w, p, msh)

    #    nV = self.nVar
    #    nMat = nV*msh.nP
    #    jac = zeros([nMat,nMat])

    #    jacL = zeros([nV,nV])
    #    jacR = zeros([nV,nV])
 
    #    wPert = w.copy()
    #    wAuxPert = wAux.copy()

    #    for iV in range(self.nVar):

    #        for iP in range(msh.nP):
    #        # perturb variables
    #            if (dw[iV,iP] > 0.):
    #               epsFD = 1./max(epsMSqrt,delta*dw[iV,iP])
    #            else:
    #               epsFD = 1./min(-epsMSqrt,delta*dw[iV,iP])

    #            wPert[iV,iP] += eps
    #            wAuxPert[iV,iP] = self.updateExtraVariables(wPert[iV,iP])

    #            iF = hstack((msh.egCloud0[iP],msh.egCloud1[iP])) # indices of flux points

    #            # get perturbed solutions at flux points 
    #            # (only return those on sides of the faces of interest.)
    #            wLP,wRP,fLP,fRP = spaceSchm.spSol2FpFlx(self, wPert, wAuxPert, 
    #                                                    msh, iF) 

    #            #wL[] = wLP? ...
 
    #            fN = zeros([self.nVar,msh.nNbrs[iP]])
    #            for ifs in self.__interfaceSolvers:
    #                fN += ifs.interfaceNormalFlux(wL,wR,fL,fR,msh.sF,msh,
    #                                              vstack([w,p]), \
    #                                              'interior')

    #            dwP = (spaceSchm.fpFlx2SpRes(self,fN,msh,iP,faceList)

    #            jacL[msh.eCloud0[iP]] \
    #                += (tile(dwP[:,0]) - dw[:,iP])*epsFD
    #            jacR[msh.eCloud1[iP]] \
    #                -= (dwP[:,1:] - dw[:,iP])*epsFD

    #            # add to petsc here. (use sciPy for now)
    #            jac[] = jacL
    #            jac[] = jacR

    #            wPert[iV,iP] = w[iV,iP] # reset variable
    #            wAuxPert[iV,iP] = wAux[iV,iP] # reset variable

    #    return jac, dw
    #---------------------------------------------------------------------------
    #def computeResidual(self,spaceSchm,msh,w):
    def computeResidual(self,spaceSchm):
    # normal-flux integral for finite volume
        
        msh = self.__msh
        w = self.__w
        p = self.__p

        wL,wR,fL,fR = spaceSchm.spSol2FpFlx(self, w, p, msh) 

        fN = zeros([self.nVar,msh.nE])
        for ifs in self.__interfaceSolvers:
            fN += ifs.interfaceNormalFlux(wL,wR,fL,fR,
                                          msh.sF,msh,
                                          vstack([w,p]), \
                                          'interior')

        dw = spaceSchm.fpFlx2SpRes(self,fN,msh)

        #print 'dw pre boundary'
        #for iV in range(self.nVar):
        #    print dw[iV,:]

        for iH in range(msh.numPatches()):
            dw[:,msh.facesOnPatch(self.__boundaryConditions[iH].patchNum())]\
         += self.__boundaryConditions[iH].faceNormalFlux()
         #+= self.__boundaryConditions[iH].faceNormalFlux(w,self,msh)

        #print 'dw post boundary'
        #for iV in range(self.nVar):
        #    print dw[iV,:]
        return dw
    #---------------------------------------------------------------------------
    #def correctBoundaryValues(self,w,msh):
    def correctBoundaryValues(self,itType='full'):
        w = self.__w
        p = self.__p
        msh = self.__msh
        for iH in range(msh.numPatches()):
            faceList = msh.facesOnPatch( \
                           self.__boundaryConditions[iH].patchNum())
            w[:,faceList], p[:,faceList] \
          = self.__boundaryConditions[iH]\
                .faceValueCorrection(itType=itType)
          #= self.__boundaryConditions[iH].faceValueCorrection(w,p,self,msh)

        #self.updateP()
        #return w
    #---------------------------------------------------------------------------
    #def updateSolutions(self, w, dw, dt, msh):
    def updateSolution(self,dw,dt):
        msh = self.__msh
        w = self.__w
        p = self.__p
        w[:,...] = w[:,...] - dw*tile(dt/msh.mass,[self.nVar,1])
        #w[:-1,...] = w[:-1,...] - dw*tile(dt/msh.mass,[self.nVar,1])

        if any(w[0,:] < 0.):
            print ' *** WARNING: negative density set to zero during update ***'
            w[0,(w[0,:]<0.)] = 0.

        self.updateP()

        #return w
   
    #---------------------------------------------------------------------------
    #def interfaceNormalFlux(self,wL,wR,fL,fR,sF,msh):

    #    sigma0 = 1  # entropy correction parameter
    #    nF = sF.shape[1]
    #    fN = zeros([self.nVar,nF])
    #    fDis = zeros([self.nVar,nF])
    #    #raw_input()
    #    normSF = sqrt(sum(sF**2,0))
    #    sF0 = sF / tile(normSF,[2,1])
    #    gma = self.__gma
    #    gm1 = self.__gma - 1
    #    
    #    for iF in range(nF):

    #        # copy face components to scalars
    #        sFF0 = sF0[0,iF]
    #        sFF1 = sF0[1,iF]

    #        # compute necessary left and right states
    #        rhoL = wL[0,iF]
    #        rhoR = wR[0,iF]
    #        rhouL = wL[1,iF]
    #        rhouR = wR[1,iF]
    #        rhovL = wL[2,iF]
    #        rhovR = wR[2,iF]
    #        pL = wL[4,iF]
    #        pR = wR[4,iF]

    #        uL = rhouL / rhoL
    #        uR = rhouR / rhoR
    #        vL = rhovR / rhoL
    #        vR = rhovR / rhoR
    #        cL = sqrt(gma*pL/rhoL)
    #        cR = sqrt(gma*pR/rhoR)

    #        # compute Roe average face quantities
    #        rhoLRt = sqrt(rhoL)
    #        rhoRRt = sqrt(rhoR)
    #        deNom = rhoLRt + rhoRRt
    #        rhoF = rhoLRt * rhoRRt;
    #        uF = (rhouL/rhoLRt + rhouR/rhoRRt) / deNom
    #        vF = (rhovL/rhoLRt + rhovR/rhoRRt) / deNom
    #        hF = ((wL[3,iF] + pL)/rhoLRt + (wR[3,iF] + pR)/rhoRRt) \
    #           / deNom
    #        velMagSq = (uF**2 + vF**2)
    #        cFSq = (gm1)*(hF - .5*(velMagSq))
    #        cF = sqrt(cFSq)

    #        velNF = uF*sFF0 + vF*sFF1
    #        velNL = uL*sFF0 + vL*sFF1
    #        velNR = uR*sFF0 + vR*sFF1

    #        # eigenvalues
    #        lamb = array([velNF,velNF,velNF+cF,velNF-cF])
    #        lambL = array([velNL,velNL,velNL+cL,velNL-cL])
    #        lambR = array([velNR,velNR,velNR+cR,velNR-cR])

    #        ## correct wave speed
    #        epsF = sigma0*maximum(zeros(self.nVar), \
    #                              lambL - lamb, lambR - lamb)
    #        ## epsF = sigma0*max(0,velNL-velNF, velNR-velNF)
    #        #
    #        #corrSpeed = [ii for ii in abs(lamb) < epsF]
    #        iCorr = arange(self.nVar)[abs(lamb)<epsF]
    #        #print lamb
    #        #if (any(corrSpeed)):
    #        #   print 'wave speed correction for face ', iF
    #        #   print msh.x[msh.e0[iF]]
    #        #   print msh.x[msh.e1[iF]]
    #        #   print uF, vF, sFF0, sFF1, velNF, cF
    #        #   print lamb
    #        #   print corrSpeed
    #        #   print lamb[iCorr]
    #        lamb[iCorr] = .5*(lamb[iCorr]**2/epsF[iCorr] + epsF[iCorr])
    #        #else:
    #        #   print 'no wave speed correction for face ', iF
    #        #print lamb

    #        ## extra smoothing when lamda is almost zero
    #        #epsSq = 1e-10
    #        #iPos = arange(self.nVar)[lamb>0]
    #        #iNeg = arange(self.nVar)[lamb<0]
    #        #print arange(self.nVar)
    #        #print iPos, iNeg
    #        #print sqrt(lamb[iPos]**2 + epsSq)
    #        #print sqrt(lamb[iNeg]**2 + epsSq)
    #        #lamb[iPos] = .5*(lamb[iPos] + sqrt(lamb[iPos]**2 + epsSq))
    #        #lamb[iNeg] = .5*(lamb[iNeg] - sqrt(lamb[iNeg]**2 + epsSq))
    #        #print 'after smoothing'
    #        #print lamb
    #        
    #        # build right and left eigenvectors
    #        bta = .5/cFSq
    #        phiSq = .5*gm1*velMagSq
    #        s0c = sFF0*cF
    #        s1c = sFF1*cF
    #        velNc = velNF*cF
    #        gm1u = gm1*uF
    #        gm1v = gm1*vF
    #        rMat = array([[1, 0, bta, bta],
    #                     [uF,  sFF1, bta*(uF+s0c), bta*(uF-s0c)],
    #                     [vF, -sFF0, bta*(vF+s1c), bta*(vF-s1c)],
    #                     [phiSq/gm1, sFF1*uF-sFF0*vF,\
    #                           bta*(hF+velNc), bta*(hF-velNc)]])
    #        lMat = array([[(1-phiSq/cFSq), gm1u/cFSq, gm1v/cFSq, (-gm1/cFSq)],
    #                     [-sFF1*uF+sFF0*vF, sFF1, -sFF0, 0.],
    #                     [phiSq-velNc,  s0c-gm1u,  s1c-gm1v, gm1],
    #                     [phiSq+velNc, -s0c-gm1u, -s1c-gm1v, gm1]])

    #        dwF = wR[:-1,iF] - wL[:-1,iF]

    #        #fNL = fL[0,:,iF] * sFF0 + fL[1,:,iF] * sFF1
    #        #fNR = fR[0,:,iF] * sFF0 + fR[1,:,iF] * sFF1 
    #        fNC = ( (fL[0,:,iF] + fR[0,:,iF]) * sFF0
    #               +(fL[1,:,iF] + fR[1,:,iF]) * sFF1) * normSF[iF]

    #        fDis[:,iF] = -dot(
    #                          dot(
    #                              rMat*tile(abs(lamb),[self.nVar,1]),
    #                              lMat
    #                             ),
    #                          dwF)\
    #                   * normSF[iF]
    #        #print ADis

    #        #fN[:,iF] = .5*((fNL + fNR)*normSF[iF] + fDis[:,iF])
    #        fN[:,iF] = .5*(fNC + fDis[:,iF])

    #    #print fN
    #    #print 'fDis'
    #    #print fDis
    #    #print 'fN'
    #    #print fN
    #    return fN
    #---------------------------------------------------------------------------
##------------------------------------------------------------------------------
## end class EulerEquations
##------------------------------------------------------------------------------

