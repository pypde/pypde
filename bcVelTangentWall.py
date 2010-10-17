import sys
from numpy import *
import bcPriority

class BCVelTangentWall(object):
    priority = bcPriority.bcPriorityDict.get('BCVelTangentWall')
    ##--------------------------------------------------------------------------
    def __init__(self,eqn='',msh='',patchNum='',refState='',
                 momentumCorrection='projection',correctWithGrad=0):
        self.__eqn = eqn
        self.__msh = msh
        self.__refState = refState
        self.__patchNum = patchNum
        self.__startFace = msh.startFaceOfPatch(patchNum)
        self.__endFace = msh.endFaceOfPatch(patchNum)
        self.__numFaces = msh.numFacesOnPatch(patchNum)
        self.__correctWithGrad = correctWithGrad
        self.__faceCorrMode = momentumCorrection
        if momentumCorrection == 'projection':
            self.faceValueCorrection = self.__projectWallVelocity
        elif momentumCorrection == 'rotation':
            self.faceValueCorrection = self.__rotateWallVelocity
        else:
            sys.exit(" *** Error: Unsupported wall face value correction mode!")
 
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
        p = eqn.p()[nodeInd]
        fN = zeros([4,self.__numFaces])
        fN[1:3,:] = tile(p,[2,1]) * msh.sB[:,nodeInd]
        return fN
    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    def __rotateWallVelocity(self,itType=None):

        nodeInd = arange(self.__startFace, self.__endFace)
        eqn = self.__eqn
        msh = self.__msh

        w = eqn.solution()
        nVar = eqn.nVar

        normSB = sqrt(sum(msh.sB[:,nodeInd]*msh.sB[:,nodeInd],0))
        sB0 = msh.sB[:,nodeInd] / tile(normSB,[2,1])
        momNorm = sqrt(sum(w[1:3,nodeInd]**2,0))
        momDotN = sum(w[1:3,nodeInd] * sB0,0)
        #w[1:3,:] -= tile(momDotN,[2,1]) * sB0
        wT = w[:,nodeInd]

        # project momentum to face tangent
        wT[1:3,:] -= tile(momDotN,[2,1]) * sB0
        momNormNew =  sqrt(sum(wT[1:3,nodeInd]**2,0))

        rotInd = arange(self.__numFaces)[momNormNew > 1e-6]
        momRatio = zeros([1,self.__numFaces])
        momRatio[0,rotInd] = momNorm[rotInd] / momNormNew[rotInd]
        wT[1:3,:] *= tile(momRatio,[2,1])

        return wT, eqn.updateP(wT)
    ##--------------------------------------------------------------------------
    #def faceValueCorrection(self,w,eqn,msh):
    #def faceValueCorrectionProjec(self):
    def __projectWallVelocity(self,itType=None):

        nodeInd = arange(self.__startFace, self.__endFace)
        eqn = self.__eqn
        msh = self.__msh

        w = eqn.solution()
        nVar = eqn.nVar
        #wT = zeros([nVar,self.__numFaces])
        #if self.__correctWithGrad:
        #    for iP2 in range(self.__numFaces):
        #        iP = nodeInd[iP2]
        #        nNb = msh.nNbr[iP]
        #        c0 = msh.clouds0[iP]
        #        c1 = msh.clouds1[iP]
        #        ec0 = msh.egClouds0[iP]
        #        ec1 = msh.egClouds1[iP]
        #        dxc = -1./3./nNb *(  sum(msh.x[c0,:],0) + sum(msh.x[c1,:],0) \
        #                           - msh.x[iP,:])
        #        #print msh.sF[:,ec0].shape
        #        #print msh.sF[:,ec1].shape
        #        #print w[:,c0].shape
        #        #print w[:,c0]
        #        #print w[:,c1].shape
        #        #print w[:,c1]
        #        ind = hstack([iP,c0,c1])
        #        gradWLoc = dot(
        #                       w[:,hstack([iP,c0,c1])], \
        #                       hstack([ reshape(msh.sB[:,iP],[2,1]),\
        #                                msh.sF[:,ec0],\
        #                               -msh.sF[:,ec1]]).transpose()\
        #                      )\
        #                 / (2*msh.mass[iP])
        #        dwLoc = sum(gradWLoc*tile(dxc,[nVar,1]),1)
        #        print dwLoc
        #        sys.stdout.flush()
        #        #raw_input()
        #        wT[:,iP2] = w[:,iP] + dwLoc

        #    #print wT.shape
        #    #print w[:,nodeInd].shape
        #    w[:,nodeInd] = wT

        normSB = sqrt(sum(msh.sB[:,nodeInd]*msh.sB[:,nodeInd],0))
        sB0 = msh.sB[:,nodeInd] / tile(normSB,[2,1])
        momDotN = sum(w[1:3,nodeInd] * sB0,0)
        #w[1:3,:] -= tile(momDotN,[2,1]) * sB0
        wT = w[:,nodeInd]

        # project momentum to face tangent
        wT[1:3,:] -= tile(momDotN,[2,1]) * sB0

        return wT, eqn.updateP(wT)
    ##--------------------------------------------------------------------------
