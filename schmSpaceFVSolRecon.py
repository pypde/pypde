from schmLocalReconConstructor import *
#from numpy import *
#from numpy.linalg import *

##------------------------------------------------------------------------------
class SchmSpaceFVSolRecon(object):
    #---------------------------------------------------------------------------
    # constructor
    #
    def __init__(self,eqn,msh,localReconDict):

        self.__schmLocalRecon \
            = constructSchmLocalRecon(eqn,msh,localReconDict)
    #---------------------------------------------------------------------------
    def localReconSchm(self):
        return self.__schmLocalRecon
    #---------------------------------------------------------------------------
    #def spSol2FpFlx(self,eqn,wDict,wAuxDict,msh):
    def spSol2FpFlx(self,eqn,w,wAux,msh):

        wPriL,wPriR = self.__schmLocalRecon.sp2Fp(w,msh)
        wAuxL,wAuxR = self.__schmLocalRecon.sp2Fp(wAux,msh)
        fL = eqn.pointFlux(wPriL,wAuxL) # needs dictionary
        fR = eqn.pointFlux(wPriR,wAuxR)
        return vstack([wPriL,wAuxL]),vstack([wPriR,wAuxR]),fL,fR 
    #---------------------------------------------------------------------------
    #def selectedSpSol2FpFlx(self,eqn,w,wAux,msh,lSP):
    ##def spSol2FpFlx(self,eqn,w,wAux,msh):
    #    
    #    f = eqn.pointFlux(w,wAux)

    #    # get list of flux points (faces) involved with points in lSP
    #    lFP = [msh.egCloud0 for iP in lSP]
    #      
    #    fL,fR = self.__schmLocalRecon.sp2SelectedFp(f,msh,iF)
    #    wL,wR = self.__schmLocalRecon.sp2SelectedFp(vstack([w,wAux]),msh,iF)
    #    return wL,wR,fL,fR 
    #---------------------------------------------------------------------------
    def fpFlx2SpRes(self,eqn,fN,msh): #,iP=None):

        dw = empty([eqn.nVar,msh.nP])
        for iV in range(eqn.nVar):
            dw[iV,:] = \
                array([sum(jj,0)
                           for jj in
                                zip([sum(fN[iV,ii]) for ii in msh.egClouds0],
                                    [-sum(fN[iV,ii]) for ii in msh.egClouds1])
                          ])

        return dw
    #---------------------------------------------------------------------------
##------------------------------------------------------------------------------
## end class SchmFluxRecons
##------------------------------------------------------------------------------
