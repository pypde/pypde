from schmLocalReconConstructor import *
#from numpy import *
#from numpy.linalg import *

##------------------------------------------------------------------------------
class SchmSpaceFVFluxRecon(object):
    #---------------------------------------------------------------------------
    # constructor
    #
    def __init__(self,msh,localReconDict):

        self.__schmLocalRecon \
            = constructSchmLocalRecon(msh,localReconDict)
    #---------------------------------------------------------------------------
    def localReconSchm(self):
        return self.__schmLocalRecon
    #---------------------------------------------------------------------------
    def spSol2FpFlx(self,eqn,w,wAux,msh,iF=None):
    #def spSol2FpFlx(self,eqn,w,wAux,msh):
        
        f = eqn.pointFlux(w,wAux)
        fL,fR = self.__schmLocalRecon.sp2Fp(f,msh,iF)
        wL,wR = self.__schmLocalRecon.sp2Fp(vstack([w,wAux]),msh,iF)
        return wL,wR,fL,fR 
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

        #if iP==None:
           
        #dw = zeros([eqn.nVar,msh.nP])
        #for e in range(msh.nE):
        #    dw[:,msh.e0[e]] += fN[:,e]  # aijs are embedded in fN\
        #    dw[:,msh.e1[e]] -= fN[:,e]

        #dw = empty([eqn.nVar,msh.nP])
        #dw = array([sum(jj,0) for jj in
        #               zip([sum(fN[:,ii]) for ii in msh.egClouds0],
        #                   [-sum(fN[:,ii]) for ii in msh.egClouds1])])

        dw = empty([eqn.nVar,msh.nP])
        for iV in range(eqn.nVar):
            dw[iV,:] = \
                array([sum(jj,0)
                           for jj in
                                zip([sum(fN[iV,ii]) for ii in msh.egClouds0],
                                    [-sum(fN[iV,ii]) for ii in msh.egClouds1])
                          ])
        # unused timing snippet
        #.......................................................................
        ##nPeat = 100
        ##t1 = timeit.Timer("dw = dw*0.\nfor e in range(msh.nE):\n    dw[:,msh.e0]] += fN[:,e]  # aijs are embedded in fN\n    dw[:,msh.e1[e]] -= fN[:,e]", \
        ##                  "from numpy import *\ndw = zeros([4,msh.nP])", \
        ##                      'gc.enable()').repeat(nPeat)
        ##print float(sum(t1)) / len(t1)
        ##sys.exit(0)
        #.......................................................................

        #else
        #    dw = empty([eqn.nVar,1])
        #    for iV in range(eqn.nVar):
        #        dw[iV,:] = \
        #            array([sum(jj,0)
        #                       for jj in
        #                            zip([sum(fN[iV,ii]) for ii in msh.egClouds0[iP]],
        #                                [-sum(fN[iV,ii]) for ii in msh.egClouds1[iP]])
        #                  ])

        return dw
    #---------------------------------------------------------------------------
##------------------------------------------------------------------------------
## end class SchmFluxRecons
##------------------------------------------------------------------------------
