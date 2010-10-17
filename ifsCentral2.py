#from numpy import *
#from numpy.linalg import *
#from pysparse import spmatrix

##------------------------------------------------------------------------------
class InterfaceSolverCentral2(object):
    #---------------------------------------------------------------------------
    # constructor
    #
    typeName = 'central2'
    def __init__(self,msh):
        # actually this is supposed to do nothing
        self.nD = 2   # 2D
    #---------------------------------------------------------------------------
    def interfaceNormalFlux (self,wL,wR,fL,fR,sF,msh,w=None,mode='interior'):
    # fL: nD x nVar x nF
        fNC =.5*( (fL[0,:,:] + fR[0,:,:]) * sF[0,:]
                 +(fL[1,:,:] + fR[1,:,:]) * sF[1,:])

        return fNC
    #---------------------------------------------------------------------------
##------------------------------------------------------------------------------
## end class LocalReconConstant
##------------------------------------------------------------------------------
