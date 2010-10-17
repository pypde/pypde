import os
from numpy import *

class ProbeForce(object):
    ##--------------------------------------------------------------------------
    def __init__(self, prefix,
                       refState='', 
                       msh='', 
                       name='',
                       patches='',  
                       lRef=1., xCR=[0.,0.],
                       writeForces=0, writeForceCoeffs=1):
                       #cleanDirectory=0):

        if writeForces == 1:
            sys.exit(' *** ERROR: Only force coefficients supported for force probes. ***')
        if writeForceCoeffs == 0:
            print(' *** WARNING: Only force coefficients supported for force probes. ***')
            print(' ***          Setting force probe to output force coefficients.   ***')
            writeForceCoeffs = 1

        self.__refState = refState
        self.__lRef = float(lRef)
        self.__xCR = array(xCR)
        #self.__startFace = msh.startFaceOfPatch(patchNum)
        #self.__endFace = msh.endFaceOfPatch(patchNum)
        #self.__numFaces = msh.numFacesOnPatch(patchNum)
        self.__forcePatchesNumList = sorted([msh.patchNum(iPatch) for iPatch in patches])
        self.__numFaces = sum(int(msh.numFacesOnPatch(iP))
                               for iP in self.__forcePatchesNumList)
        # create or clean directory
        self.__probeDir = ''.join([prefix,'/',name])
        if not(os.path.isdir(self.__probeDir)):
            #if cleanDirectory && (prefix != None):
                 #print probeDir
	    #     #os.system(''.join(['rm ',probeDir,'/*']))
        #else: 
            os.mkdir(self.__probeDir)
        self.__filename = ''.join([self.__probeDir,'/',name,'.forceCoeffs'])
        fileProbe = open(self.__filename,'w')
        fileProbe.write('#it t cL cD cM\n')
        fileProbe.close()
       
    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    def calcPressureForces(self,w,msh):
        

        alphaRef = self.__refState.alpha()*pi/180.

        #nodeInd = arange(self.__startFace, self.__endFace)
        nodeInd = hstack([arange(msh.startFaceOfPatch(iP),
                                      msh.endFaceOfPatch(iP)) 
                               for iP in self.__forcePatchesNumList])

        p = w[4,nodeInd]
        rhoRef = self.__refState.rho()
        velRef = self.__refState.velMag()
        qInfF = .5*rhoRef*(velRef**2) * self.__lRef
        cA = cos(alphaRef)
        sA = sin(alphaRef)

        xR = msh.x[nodeInd,:] - tile(self.__xCR,[self.__numFaces,1])
        xR[:,1] = -xR[:,1]
        xR = fliplr(xR)  # xR is [y -x] for moment calculation

        pF = msh.sB[:,nodeInd] * tile(p,[2,1])
        forceVec = sum(pF,1)
        lift = sum(forceVec*array([-sA,cA]))
        drag = sum(forceVec*array([cA,sA]))
        #circ = lift/(rhoRef * velRef)
 
        #print pF.shape
        #print xR.shape
        #print msh.sB.shape

        # only 1 component for 2-D moment
        mmnt = sum(pF*xR.transpose())

        cL = lift / qInfF
        cD = drag / qInfF
        cM = mmnt / (qInfF*self.__lRef)

        return cL, cD, cM
    ##--------------------------------------------------------------------------
    def writeProbeValues(self,it,w,msh):
        p = w[4,:]	# HACK!
        cL, cD, cM = self.calcPressureForces(w,msh)

        fileProbe = open(self.__filename,'a')
        fileProbe.write('%d %24.12e %24.12e %24.12e\n' % (it,cL,cD,cM))
        fileProbe.close()

    ##--------------------------------------------------------------------------
    def writePressure(self,it,w,msh):
        nodeInd = hstack([arange(msh.startFaceOfPatch(iP),
                                 msh.endFaceOfPatch(iP)) 
                          for iP in self.__forcePatchesNumList])

        
        #p = reshape(w[4,nodeInd],[self.__numFaces,1])
        #print p
        p = w[4,nodeInd]

        pRef = self.__refState.p()
        rhoRef = self.__refState.rho()
        velRef = self.__refState.velMag()
        qInf = .5*rhoRef*(velRef*velRef)

        cp = (p-tile(pRef,p.shape))/qInf

        fnameCP = ''.join([self.__probeDir,'/','cp_',str(it),'.txt'])
        fileCP = open(fnameCP,'w')
        fileCP.write('# x y cp\n')
        xB = msh.x[nodeInd,:]
        #(hstack([ msh.x[nodeInd,:], \
        #         ((p-tile(pRef,p.shape))/qInf)])\
        #).tofile(fileCP,format='%24.12e',sep=' ')
        #fileCP.write('\n')
        for iF in range(self.__numFaces):
            fileCP.write('%24.12e %24.12e %24.12e\n' %(xB[iF,0], xB[iF,1], cp[iF]))
        fileCP.close()
