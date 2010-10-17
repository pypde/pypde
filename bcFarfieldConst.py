from numpy import *
from ifsConstructor import *
import bcPriority

class BCFarfieldConst(object):
    priority = bcPriority.bcPriorityDict.get('BCFarfieldConst')
    ##--------------------------------------------------------------------------
    def __init__(self,eqn=None,msh=None,patchNum=None,refState=None,ifsDict=None):
        #self.__refState = refState
        self.__patchNum = patchNum
        self.__startFace = msh.startFaceOfPatch(patchNum)
        self.__endFace = msh.endFaceOfPatch(patchNum)
        self.__numFaces = msh.numFacesOnPatch(patchNum)
        self.__interfaceSolvers = constructInterfaceSolver(eqn,msh,ifsDict)
        rho0 = refState.rho()
	u0 = refState.u()
	v0 = refState.v()
        p0 = refState.p()
        e0 = rho0*refState.h() - p0
        self.__w0 = tile([rho0,rho0*u0, rho0*v0,rho0*e0],[self.__numFaces,1]).transpose()
        self.__p0 = tile(p0,[1,self.__numFaces])
        self.__f0 = eqn.pointFlux(self.__w0, self.__p0)
        self.__eqn = eqn
        self.__msh = msh
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

    ##--------------------------------------------------------------------------
    #def faceValueCorrection(self,w,p,eqn,msh):
    def faceValueCorrection(self,itType=None):
        #return self.__w0, self.__p

        nodeInd = arange(self.__startFace, self.__endFace)
        eqn = self.__eqn
        msh = self.__msh
        wCons = eqn.solution()[:,nodeInd]
        p = eqn.p()[nodeInd]

        return wCons, p  # i.e. do nothing!
    ##--------------------------------------------------------------------------
