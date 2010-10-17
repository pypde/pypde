from numpy import *
#from numpy.linalg import *
#from pysparse import spmatrix

##------------------------------------------------------------------------------
class InterfaceSolverRoe(object):
    #---------------------------------------------------------------------------
    # constructor
    #
    typeName = 'roe'
    def __init__(self,msh=None, sigma=1.,gma=1.4):
        self.__sigma0 = sigma
        self.__gma = gma
    #---------------------------------------------------------------------------
    def interfaceNormalFlux (self,wL,wR,fL,fR,sF,msh=None,w=None,mode='interior'):
        """ Computes Roe dissipative flux (area weighted) for Euler equations.

        *** CONTRIBUTION TO DISSIPATION ONLY!  
            CENTRAL PART MUST BE SPECIFIED SEPARATELY! ***

        Asssumptions:
        - CONSTANT PROPERTY GAS

        Input arguments:
        - wL:  variables values (including auxiliary ones) on left of face
               (nVar x nF)
        - wR:  variables values (including auxiliary ones) on right of face
               (nVar x nF)
        - fL:  flux values (only for conservative variables) on left of face
               (nD x nVar x nF)
        - fR:  flux values (only for conservative variables) on right of face
               (nD x nVar x nF)
        - sF:  array of face vectors (nD x nF)
        - w :  pointwise variable values (including auxiliary ones) (nVar x nP)
               ***NOT USED***
        - mode :  flag to select code for interior or boundary calculations

        Output arguments:
        - fN:  Roe dissipative flux (area weighted) on each face

        Notes:
        - conservative variables (order expected in w):  rho, rhou, rhov, rhoe
        - auxiliary variables:  p
        - number of spatial dimensions: nD
        - number of faces: nF
        - number of points: nP
        - number of equation (conservative) variables: nVar

        """


        sigma0 = self.__sigma0  # entropy correction parameter
        nF = wL.shape[1]
        #nF = msh.nE
        #print msh.nE
        #print fL.shape
        nVar = fL.shape[1]

        fN = zeros([nVar,nF])
        normSF = sqrt(sum(sF**2,0))
        sF0 = sF / tile(normSF,[2,1])
        gma = self.__gma
        gm1 = self.__gma - 1
        
        for iF in range(nF):

            # copy face components to scalars
            sFF0 = sF0[0,iF]
            sFF1 = sF0[1,iF]

            # compute necessary left and right states
            rhoL = wL[0,iF]
            rhoR = wR[0,iF]
            rhouL = wL[1,iF]
            rhouR = wR[1,iF]
            rhovL = wL[2,iF]
            rhovR = wR[2,iF]
            pL = wL[4,iF]
            pR = wR[4,iF]

            uL = rhouL / rhoL
            uR = rhouR / rhoR
            vL = rhovR / rhoL
            vR = rhovR / rhoR
            cL = sqrt(gma*pL/rhoL)
            cR = sqrt(gma*pR/rhoR)

            # compute Roe average face quantities
            rhoLRt = sqrt(rhoL)
            rhoRRt = sqrt(rhoR)
            deNom = rhoLRt + rhoRRt
            rhoF = rhoLRt * rhoRRt;
            uF = (rhouL/rhoLRt + rhouR/rhoRRt) / deNom
            vF = (rhovL/rhoLRt + rhovR/rhoRRt) / deNom
            hF = ((wL[3,iF] + pL)/rhoLRt + (wR[3,iF] + pR)/rhoRRt) \
               / deNom
            velMagSq = (uF**2 + vF**2)
            cFSq = (gm1)*(hF - .5*(velMagSq))
            cF = sqrt(cFSq)

            velNF = uF*sFF0 + vF*sFF1
            velNL = uL*sFF0 + vL*sFF1
            velNR = uR*sFF0 + vR*sFF1

            # eigenvalues
            lamb = array([velNF,velNF,velNF+cF,velNF-cF])
            lambL = array([velNL,velNL,velNL+cL,velNL-cL])
            lambR = array([velNR,velNR,velNR+cR,velNR-cR])

            ## correct wave speed
	    epsF = sigma0*maximum(zeros(nVar), \
                                  lambL - lamb, lambR - lamb)
	    ## epsF = sigma0*max(0,velNL-velNF, velNR-velNF)
            #
            #corrSpeed = [ii for ii in abs(lamb) < epsF]
            iCorr = arange(nVar)[abs(lamb)<epsF]
            #print lamb
            #if (any(corrSpeed)):
            #   print 'wave speed correction for face ', iF
            #   print msh.x[msh.e0[iF]]
            #   print msh.x[msh.e1[iF]]
            #   print uF, vF, sFF0, sFF1, velNF, cF
            #   print lamb
            #   print corrSpeed
            #   print lamb[iCorr]
            lamb[iCorr] = .5*(lamb[iCorr]**2/epsF[iCorr] + epsF[iCorr])
            #else:
            #   print 'no wave speed correction for face ', iF
            #print lamb

            ## extra smoothing when lamda is almost zero
            #epsSq = 1e-10
            #iPos = arange(nVar)[lamb>0]
            #iNeg = arange(nVar)[lamb<0]
            #print arange(nVar)
            #print iPos, iNeg
            #print sqrt(lamb[iPos]**2 + epsSq)
            #print sqrt(lamb[iNeg]**2 + epsSq)
            #lamb[iPos] = .5*(lamb[iPos] + sqrt(lamb[iPos]**2 + epsSq))
            #lamb[iNeg] = .5*(lamb[iNeg] - sqrt(lamb[iNeg]**2 + epsSq))
            #print 'after smoothing'
            #print lamb
            
            # build right and left eigenvectors
            bta = .5/cFSq
            phiSq = .5*gm1*velMagSq
            s0c = sFF0*cF
            s1c = sFF1*cF
            velNc = velNF*cF
            gm1u = gm1*uF
            gm1v = gm1*vF
            rMat = array([[1, 0, bta, bta],
                         [uF,  sFF1, bta*(uF+s0c), bta*(uF-s0c)],
                         [vF, -sFF0, bta*(vF+s1c), bta*(vF-s1c)],
                         [phiSq/gm1, sFF1*uF-sFF0*vF,\
                               bta*(hF+velNc), bta*(hF-velNc)]])
            lMat = array([[(1-phiSq/cFSq), gm1u/cFSq, gm1v/cFSq, (-gm1/cFSq)],
                         [-sFF1*uF+sFF0*vF, sFF1, -sFF0, 0.],
                         [phiSq-velNc,  s0c-gm1u,  s1c-gm1v, gm1],
                         [phiSq+velNc, -s0c-gm1u, -s1c-gm1v, gm1]])

            dwF = wR[:-1,iF] - wL[:-1,iF]

            fN[:,iF] = -dot(
                              dot(
                                  rMat*tile(abs(lamb),[nVar,1]),
                                  lMat
                                 ),
                              dwF)\
                       * normSF[iF] *.5

        return fN
##------------------------------------------------------------------------------
## end class InterfaceSolverRoe
##------------------------------------------------------------------------------
