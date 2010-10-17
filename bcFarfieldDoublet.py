from numpy import *
from ifsConstructor import *
import bcPriority

class BCFarfieldDoublet(object):
    priority = bcPriority.bcPriorityDict.get('BCFarfieldDoublet')
    ##--------------------------------------------------------------------------
    def __init__(self,eqn='',msh='',patchNum='',refState='',xCR=[0.,0.],
                 kappa = 0., rEquiv = 1., ifsDict=None):
        self.__eqn = eqn
        self.__msh = msh
        #self.__refState = refState
        self.__patchNum = patchNum
        self.__startFace = msh.startFaceOfPatch(patchNum)
        self.__endFace = msh.endFaceOfPatch(patchNum)
        self.__numFaces = msh.numFacesOnPatch(patchNum)
        self.__interfaceSolvers = constructInterfaceSolver(eqn,msh,ifsDict)

        #rho0 = refState.rho()
	u0 = refState.u()
	v0 = refState.v()
        #p0 = refState.p()
        #e0 = rho0*refState.h() - p0
        h0 = refState.h()
        s0 = refState.s()
        gma = refState.gma()
        velRef = refState.velMag()
        mach = refState.mach()
        bta = sqrt(1.-mach*mach)

        #xCR = array(xCR.asList())
        xCR = array(xCR)
        #print xCR

        if kappa == 0:
            kappa = 2*pi*rEquiv*rEquiv*velRef

        nodeInd = arange(self.__startFace, self.__endFace)
        xR = msh.x[nodeInd,:] - tile(xCR,[self.__numFaces,1])
        theta = arctan2(xR[:,1],xR[:,0])
        cT = cos(theta)
        sT = sin(theta)

        kOn2PiRSq = (kappa/(2.*pi))/sum(xR*xR,1) 
        # compute doublet velocity
        uF = u0 + kOn2PiRSq/bta*(sT*sT-cT*cT)
        vF = v0 - kOn2PiRSq*2*sT*cT

        halfqF = .5*(uF*uF + vF*vF)
        cF2 = (gma-1.) * (h0 - halfqF)
        rhoF = ((s0/gma)*cF2)**(1/(gma-1.)) 
        pF = (gma*rhoF)/cF2
        rhoeF = pF/(gma-1) + rhoF*halfqF

        self.__w0 = vstack((rhoF,rhoF*uF, rhoF*vF,rhoeF))
        self.__p0 = pF
        self.__f0 = eqn.pointFlux(self.__w0,pF)
        self.__s0 = s0
        self.__h0 = h0
        self.__gma = gma
    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    def patchNum(self):
        return self.__patchNum
    ##--------------------------------------------------------------------------
    #def faceNormalFlux(self,w,eqn,msh):
    def faceNormalFlux(self):
        nodeInd = arange(self.__startFace, self.__endFace)
        eqn = self.__eqn
        msh = self.__msh
        wCons = eqn.solution()[:,nodeInd]
        p = eqn.p()[nodeInd]
        w = vstack([wCons,p])
        #fN = eqn.interfaceNormalFlux(w,vstack([self.__w0,self.__p0]),
        #                             eqn.pointFlux(wCons,p),self.__f0,
        #                             msh.sB[:,nodeInd],msh)

        fN = zeros([eqn.nVar,self.__numFaces])
        for ifs in self.__interfaceSolvers:
            fN += ifs.interfaceNormalFlux(w,vstack([self.__w0,self.__p0]),
                                          eqn.pointFlux(wCons,p),self.__f0,
                                          msh.sB[:,nodeInd],msh,
                                          eqn.allVariables(),'boundary')
        return fN
    ##--------------------------------------------------------------------------
    #def faceValueCorrection(self,w,eqn,msh):
    def faceValueCorrection(self,itType=None):

        #return self.__w0

	nodeInd = arange(self.__startFace, self.__endFace)
        eqn = self.__eqn
        msh = self.__msh
        wCons = eqn.solution()[:,nodeInd]
        p = eqn.p()[nodeInd]

        return wCons, p

        #gma = self.__gma
        #gm1 = gma-1.
        #gmg = gm1/gma

        #nodeInd = arange(self.__startFace, self.__endFace)
        #nB = self.__numFaces
        #sB = msh.sB[:,nodeInd]
        #sB /= sqrt(sum(sB**2,0))
        #tB = vstack((-sB[1,:],sB[0,:]))
        #
        ##print tB.shape
        ##print self.__w0.shape
        #rhoI = w[0,nodeInd]
        #uI = w[1:3,nodeInd]/rhoI
        #pI = w[4,nodeInd]
        #cI = sqrt(gma*pI/rhoI)

        #u0 = self.__w0[1:3,:]/tile(self.__w0[0,:],[2,1])
        #c0 = sqrt(gma*self.__w0[4,:]/self.__w0[0,:])

        #uT = sum(u0*tB,0)
        #s = self.__s0*ones(nB)

        ## extrapolate NORMAL Riemann Invariant
        ##......................................................................
        #riI = sum(uI*sB,0) + 2./gm1*cI
        #ri0 = sum(u0*sB,0) - 2./gm1*c0

        #uN = .5*(riI + ri0)
        #ind = arange(nB)[uN > 0.]
        #s[ind] = (rhoI[ind]**gma)/pI[ind]
        #uT[ind] = sum(uI[:,ind]*tB[:,ind],0)

        #uB = tile(uN,[2,1])*sB + tile(uT,[2,1])*tB
        #halfVelMagSq = .5*(sum(uB**2,0))
        #pOnRho = gmg*(self.__h0 - halfVelMagSq)

        #rhoB = (s*pOnRho)**(1./gm1)
        #rhouB = tile(rhoB,[2,1])*uB
        #pB = rhoB*pOnRho
        #rhoeB = pB/(gma-1.) + rhoB*halfVelMagSq
        #......................................................................

        ## extrapolate UPSTREAM Riemann Invariant
        ##......................................................................
        #uN = sum(uI*sB,0)
        #ri = uN + 2./gm1*cI
        #t = ones(nB)
        #ind = arange(nB)[uN > 0.]
        #s[ind] = (rhoI[ind]**gma)/pI[ind]
        #uT[ind] = sum(uI[:,ind]*tB[:,ind],0)
        #ri[ind] = sum(u0[:,ind]*sB[:,ind],0) - 2./gm1*c0[ind]
        #t[ind] = -1.

        #b = t*ri/gm1
        #gmf = (2./gm1 + 1.)/gm1
        ##c = .5*(ri**2. + uT**2) - self.__h0
        #c = (sqrt(b**2. - gmf*(.5*(ri**2.+uT**2)-self.__h0)) + b) / gmf

        #uN = ri - (2./gm1)*(t*c)

        #uB = tile(uN,[2,1])*sB + tile(uT,[2,1])*tB
        #halfVelMagSq = .5*(sum(uB**2,0))

        #pOnRho = c**2/gma

        #rhoB = (s*pOnRho)**(1./gm1)
        #rhouB = tile(rhoB,[2,1])*uB
        #pB = rhoB*pOnRho
        #rhoeB = rhoB*self.__h0 - pB
        #......................................................................
        #return vstack((rhoB,rhouB,rhoeB,pB))
    ##--------------------------------------------------------------------------
