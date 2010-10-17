from numpy import *
from ifsConstructor import *
import bcPriority

class BCClFarfieldCirc(object):
    priority = bcPriority.bcPriorityDict.get('BCClFarfieldCirc')
    ##--------------------------------------------------------------------------
    def __init__(self,eqn=None,msh=None,patchNum=None,refState=None,
                 patches=None,
                 xCR=[0.,0.], lRef=1., cLTarget=0.,alphaMax=10.,alphaMin=-10.,
                 refRes=1e-3,updateInterval=200, relaxation=0.01,
                 ifsDict=None):
        self.__refState = refState
        self.__xCR = array(xCR)
        self.__lRef = lRef
        self.__cLTarget = cLTarget
        self.__ralphaMax = alphaMax*pi/180.
        self.__ralphaMin = alphaMin*pi/180.
        self.__refRes = refRes
        self.__updateInterval = updateInterval
        self.__relax = relaxation
        self.__counter = 0
        self.__switch = 0
        self.__patchNum = patchNum
        self.__startFace = msh.startFaceOfPatch(patchNum)
        self.__endFace = msh.endFaceOfPatch(patchNum)
        self.__numFaces = msh.numFacesOnPatch(patchNum)
        self.__forcePatchesNumList = sorted([msh.patchNum(iPatch) for iPatch in patches])
        #print self.__forcePatchesNumList
        self.__interfaceSolvers = constructInterfaceSolver(eqn,msh,ifsDict)

        self.__eqn = eqn
        self.__msh = msh

        print self.__refRes
    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    def patchNum(self):
        return self.__patchNum
    ##--------------------------------------------------------------------------
    #def faceNormalFlux(self,w,eqn,msh):
    def faceNormalFlux(self):
        

        #u0 = self.__refState.u()
        #v0 = self.__refState.v()
        #h0 = self.__refState.h()
        #s0 = self.__refState.s()
        #gma = self.__refState.gma()
        #machRef = self.__refState.mach()
        #alphaRef = self.__refState.alpha()*pi/180.

        eqn = self.__eqn
        msh = self.__msh

        nodeInd = arange(self.__startFace, self.__endFace)
        #forceNodeInd = hstack([arange(msh.startFaceOfPatch(iP),
        #                              msh.endFaceOfPatch(iP)) 
        #                       for iP in self.__forcePatchesNumList])

        wCons = eqn.solution()[:,nodeInd]
        p = eqn.p()[nodeInd]
        #pForce = eqn.p()[forceNodeInd]
        w = vstack([wCons,p])

        #rhoRef = self.__refState.rho()
        #velRef = self.__refState.velMag()

        #cA = cos(alphaRef)
        #sA = sin(alphaRef)

        #forceVec = sum(msh.sB[:,forceNodeInd] * tile(pForce,[2,1]),1)
        #lift = sum(forceVec*array([-sA,cA]),0)
        #circ = lift/(rhoRef * velRef)

        #xR = msh.x[nodeInd,:] - tile(self.__xCR,[self.__numFaces,1])
        #r = sum(xR**2,1)
        #
        #theta = arctan2(xR[:,1],xR[:,0])
        #cT = cos(theta)
        #sT = sin(theta)

        #dVel = circ / (r*(1. - machRef*(sT*cA - cT*sA)**2))
        #uF = u0 + dVel*sT
        #vF = v0 + dVel*cT
        #halfqF = .5*(uF*uF + vF*vF)
        #cF2 = (gma-1) * (h0 - halfqF)
        #rhoF = ((s0/gma)*cF2)**(1/(gma-1))
        #pF = (gma*rhoF)/cF2
        #rhoeF = pF/(gma-1) + rhoF*halfqF

        #wF = vstack((rhoF,rhoF*uF,rhoF*vF,rhoeF))
        #
        #fF = eqn.pointFlux(wF,pF)

        fN = zeros([eqn.nVar,self.__numFaces])
        #for ifs in self.__interfaceSolvers:
        #    fN += ifs.interfaceNormalFlux(w,vstack((wF,pF)),
        #                                  eqn.pointFlux(wCons,p),fF,
        #                                  msh.sB[:,nodeInd],msh,
        #                                  eqn.allVariables(),'boundary')

        for ifs in self.__interfaceSolvers:
            fN += ifs.interfaceNormalFlux(w,w,
                                          eqn.pointFlux(wCons,p),
                                          eqn.pointFlux(wCons,p),
                                          msh.sB[:,nodeInd],msh,
                                          eqn.allVariables(),'boundary')

        return fN
    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    #def faceValueCorrection(self,w,p,eqn,msh):
    def faceValueCorrection(self,itType='full'):
        #return self.__w0, self.__p
        
        nodeInd = arange(self.__startFace, self.__endFace)
        eqn = self.__eqn
        msh = self.__msh
        wCons = eqn.solution()[:,nodeInd]
        p = eqn.p()[nodeInd]

        #return wCons, p  # i.e. do nothing!

        u0 = self.__refState.u()
        v0 = self.__refState.v()
        h0 = self.__refState.h()
        s0 = self.__refState.s()
        gma = self.__refState.gma()
        machRef = self.__refState.mach()

        forceNodeInd = hstack([arange(msh.startFaceOfPatch(iP),
                                      msh.endFaceOfPatch(iP)) 
                               for iP in self.__forcePatchesNumList])

        pForce = eqn.p()[forceNodeInd]

        rhoRef = self.__refState.rho()
        velRef = self.__refState.velMag()
        alphaRef = self.__refState.alpha()*pi/180.
        cA = cos(alphaRef)
        sA = sin(alphaRef)

        forceVec = sum(msh.sB[:,forceNodeInd] * tile(pForce,[2,1]),1)
        lift = sum(forceVec*array([-sA,cA]),0)
        circ = lift/(rhoRef * velRef)

        if itType=='full':
         
            # check residual   
            print eqn.maxResiduals().max()
            if ((eqn.maxResiduals()).max() < self.__refRes):
                if (self.__switch == 0):
                    self.__switch = 1
                    self.__counter = 0
            else:
                self.__switch = 0

            print self.__switch
            print self.__counter
            print self.__updateInterval
            if (self.__switch == 1):
                if (self.__counter % self.__updateInterval == 0):
                    # update AoA due to target Cl
                    qInfF = .5*rhoRef*(velRef**2) * self.__lRef
                    cL = lift / qInfF
                    dcl = self.__cLTarget - cL
                    #self.__refState.setAlpha(alphaRef + relax*dcl*90./(pi*pi))
                    self.__refState.setAlphaRadians(
                                       max(min(self.__ralphaMax,
                                           alphaRef + self.__relax*dcl/(2.*pi)),
                                           self.__ralphaMin)
                                    )
                    print 'AoA is', self.__refState.alpha()
                self.__counter += 1

        alphaRef = self.__refState.alpha()*pi/180.
        cA = cos(alphaRef)
        sA = sin(alphaRef)


        xR = msh.x[nodeInd,:] - tile(self.__xCR,[self.__numFaces,1])
        r = sum(xR**2,1)
        
        theta = arctan2(xR[:,1],xR[:,0])
        cT = cos(theta)
        sT = sin(theta)

        dVel = circ / (r*(1. - machRef*(sT*cA - cT*sA)**2))
        uF = u0 + dVel*sT
        vF = v0 + dVel*cT
        halfqF = .5*(uF*uF + vF*vF)
        cF2 = (gma-1) * (h0 - halfqF)
        rhoF = ((s0/gma)*cF2)**(1/(gma-1))
        pF = (gma*rhoF)/cF2
        rhoeF = pF/(gma-1) + rhoF*halfqF

        wF = vstack((rhoF,rhoF*uF,rhoF*vF,rhoeF))
        #fF = eqn.pointFlux(wF,pF)

        return wF, pF
    ##--------------------------------------------------------------------------

