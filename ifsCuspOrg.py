import sys

from decimal import Decimal, getcontext

from numpy import *
#from numpy.linalg import *
#from pysparse import spmatrix

##------------------------------------------------------------------------------
class InterfaceSolverCusp(object):
    #---------------------------------------------------------------------------
    # constructor
    #
    def __init__(self,msh,refState,vis1=1.,vis2=1.,cLim=1e-5,dLim =.5):
        self.__vis1 = vis1
        self.__vis2 = vis2
        self.__gma = refState.gma()
        bLim = vis2/(msh.nP**.75)
        #print bLim
        rho0 = refState.rho()
        c0 = refState.c()
        h0 = refState.h()
        self.__wLim = bLim*rho0* array([1., c0, c0, h0])
        self.__cLim = cLim
        self.__dLim = dLim
        #self.__eqn = eqn
        #self.__msh = msh
    #---------------------------------------------------------------------------
    def interfaceNormalFlux (self,wL,wR,fL,fR,sF,msh,w=None,mode='interior',iF=None): 
        """ Computes CUSP dissipative flux (area weighted) for Euler equations.

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
        #- sF:  array of face vectors (nD x nF)
        - w :  pointwise variable values (including auxiliary ones) (nVar x nP)
        - mode :  flag to select code for interior or boundary calculations
        - iF : array of flux point indices to get fluxes for

        Output arguments:
        - fN:  CUSP dissipative flux (area weighted) on each face

        Notes:
        - conservative variables (order expected in w):  rho, rhou, rhov, rhoe
        - auxiliary variables:  p
        - number of spatial dimensions: nD
        - number of faces: nF
        - number of points: nP
        - number of equation (conservative) variables: nVar

        """
        # *** want to take the solutions w as an argument
        # because want to reuse this routine for both interior and boundary
         
        nVar = fL.shape[1]
        if mode == 'boundary':
            fN = zeros([nVar,fL.shape[2]])
        else:
            gma = self.__gma
            gm1 = gma - 1.
            gmg = gm1/gma

            if iF == None:
                nB = msh.nB
                ec0 = msh.egClouds0
                ec1 = msh.egClouds1
                e0 = msh.e0
                e1 = msh.e1
            else:
                nB = msh.nB
                ec0 = msh.egClouds0
                ec1 = msh.egClouds1
                e0 = msh.e0[iF]
                e1 = msh.e1[iF]
               

            dwFL, dwFR = self.__getGrad()

            #sF = msh.sF
            sFe0 = sF/msh.mass[e0]
            sFe1 = sF/msh.mass[e1]


            #w = eqn.solution()
            wLC = w[:-1,e0].copy()
            pL = w[4,e0].copy()
            wRC = w[:-1,e1].copy()
            pR = w[4,e1].copy()
            
            wLC[3,:] += pL  # wLC[3,:] is enthalpy
            wRC[3,:] += pR  # wRC[3,:] is enthalpy

            #for ii in range(nVar):
            #    print wL[ii,:] - wLC[ii,:]

            # for FD jacobian, get gradient involving point with modified solution
            wF = .5*(wLC + wRC)
            wFNx0 = wF*tile(sFe0[0,:],[nVar,1])
            wFNy0 = wF*tile(sFe0[1,:],[nVar,1])

            wFNx1 = wF*tile(sFe1[0,:],[nVar,1])
            wFNy1 = wF*tile(sFe1[1,:],[nVar,1])
            #print wFNx.shape

            #for iV in range(nVar):
	    #    print array([sum(jj,0)
            #                  for jj in
            #                      zip([ sum(wFNx[iV,ii]) for ii in ec0],
            #                          [-sum(wFNx[iV,ii]) for ii in ec1])
            #                      ]).shape
            wx2 = zeros([nVar,msh.nP])
            for e in range(msh.nE):
                wx2[:,e0[e]] += w[:-1,e1[e]]*sFe0[0,e]
                wx2[3,e0[e]] += w[4,e1[e]]*sFe0[0,e]
                wx2[:,e1[e]] -= w[:-1,e0[e]]*sFe1[0,e]
                wx2[3,e1[e]] -= w[4,e0[e]]*sFe1[0,e]
            wy2 = zeros([nVar,msh.nP])
            for e in range(msh.nE):
                wy2[:,e0[e]] += w[:-1,e1[e]]*sFe0[1,e]
                wy2[3,e0[e]] += w[4,e1[e]]*sFe0[1,e]
                wy2[:,e1[e]] -= w[:-1,e0[e]]*sFe1[1,e]
                wy2[3,e1[e]] -= w[4,e0[e]]*sFe1[1,e]

            #wxC = empty([nVar,msh.nP])
            #wyC = empty([nVar,msh.nP])
            #for iV in range(nVar):
	    #    wxC[iV,:] = array([sum(jj,0)
            #                  for jj in
            #                      zip([ sum(wFNx0[iV,ii]) for ii in ec0],
            #                          [-sum(wFNx1[iV,ii]) for ii in ec1])
            #                  ])

	    #    wyC[iV,:] = array([sum(jj,0)
            #                  for jj in
            #                      zip([ sum(wFNy1[iV,ii]) for ii in ec0],
            #                          [-sum(wFNy1[iV,ii]) for ii in ec1])
            #                  ])

	    #cx  = array([sum(jj,0)
            #                for jj in
            #                    zip([ sum(sF[0,ii]) for ii in ec0],
            #                        [-sum(sF[0,ii]) for ii in ec1])
            #                ])
            #cx[:msh.nB] += msh.sB[0,:]
            #print 'cx'
            #print cx/msh.mass

            #for iV in range(nVar):
            #    print iV
            #    print 'wLC'
            #    print wLC[iV,:]
            #    print 'wRC'
            #    print wRC[iV,:]
            #    print 'w'
            #    print w[iV,:]
            #    print 'wx preCorr'
            #    print wx[iV,:]

            # correct boundary values (need to get original w somehow
            sB = msh.sB/tile(msh.mass[:nB],[2,1])
            wx2[:,:nB] += w[:-1,:nB]*tile(sB[0,:],[nVar,1])
            wy2[:,:nB] += w[:-1,:nB]*tile(sB[1,:],[nVar,1])
            wx2[3,:nB] += w[4,:nB]*sB[0,:]
            wy2[3,:nB] += w[4,:nB]*sB[1,:]

            #wxC[:,:nB] += w[:-1,:nB]*tile(sB[0,:],[nVar,1])
            #wyC[:,:nB] += w[:-1,:nB]*tile(sB[1,:],[nVar,1])
            #wxC[3,:nB] += w[4,:nB]*sB[0,:]
            #wyC[3,:nB] += w[4,:nB]*sB[1,:]

            wx = wx2*.5 # (b/c meshless aij is sF/2)
            wy = wy2*.5

            #print wx2.shape
            #print msh.mass.shape
            #wx /= msh.mass
            #wx3 = wx2/msh.mass
            #wx2 /= tile(msh.mass,[nVar,1])
            #wy /= msh.mass

            #for iV in range(nVar):
            #    print iV
            #    print 'wxC'
            #    print wxC[iV,:]
            #    print 'wx2'
            #    print wx2[iV,:]
            #    #print 'wx3'
            #    #print wx3[iV,:]
            #    #print 'wx3-wx2'
            #    #print wx3[iV,:]-wx2[iV,:]
            ##print 'wy'
            ##for iV in range(nVar):
            ##    print iV
            ##    print wy[iV,:]
            #sys.exit(0)

            dx = msh.x[e1,:] - msh.x[e0,:]
            dwFL =   wx[:,e0]*(tile(dx[:,0],[nVar,1])) \
                   + wy[:,e0]*(tile(dx[:,1],[nVar,1]))
            dwFR =   wx[:,e1]*(tile(dx[:,0],[nVar,1])) \
                   + wy[:,e1]*(tile(dx[:,1],[nVar,1]))
            #dwFL = .5*(  wx[:,e0]*(tile(dx[:,0],[nVar,1])) \
            #           + wy[:,e0]*(tile(dx[:,1],[nVar,1])))
            #dwFR = .5*(  wx[:,e1]*(tile(dx[:,0],[nVar,1])) \
            #           + wy[:,e1]*(tile(dx[:,1],[nVar,1])))
            #for iV in range(nVar):
            #    for ii in range(msh.nE):
            #        print dwFL[iV,ii], dwFR[iV,ii]

            # SLIP
            wLimA = tile(self.__wLim,[msh.nE,1]).transpose()
            aShape = [nVar,msh.nE]
            a = reshape(maximum(abs(dwFL) + abs(dwFR),wLimA),[wLimA.size])
            b = reshape(abs(dwFR - dwFL),[wLimA.size])
            
            r = ones(wLimA.size)
            #print a.size
            #print r.size
            #print aShape
            #print a.shape
            #print a>0

            ind = arange(r.size)[a > 1e-6]
            #ind = arange(r.size)[a > 0.]
            #r[ind] = 1 - (b[ind]/a[ind])**3.
            r[ind] = 1 - (b[ind]/a[ind])**2.   # r -> 0, dissipation -> 1st order.

            #print arange(r.size)[r !=1.]
            #nn1 = 0
            #nn2 = 0 
            #for ii,val in enumerate(r):
            #    if val != 1.:
            #        print ii, val
            #        nn1 += 1
            #    else:
            #        print ii, val
            #        nn2 += 1
            #print 'nn1:', nn1
            #print 'nn2:', nn2
            #sys.exit(0)
            #ab = .25*reshape(r**5.,aShape)*(dwFL + dwFR)
            ab = .25*reshape(r,aShape)*(dwFL + dwFR)  # limited average of dw's
 
            wRC -= ab
            wLC += ab

            #for iV in range(nVar):
            #    for ii in range(msh.nE):
            #        print iV, wLC[iV,ii], wRC[iV,ii]
            #sys.exit(0)

            # CUSP
            # compute necessary left and right states
            rhoL = wLC[0,:]
            rhoR = wRC[0,:]
            rhouL = wLC[1,:]
            rhouR = wRC[1,:]
            rhovL = wLC[2,:]
            rhovR = wRC[2,:]
            rhohL = wLC[3,:]
            rhohR = wRC[3,:]

            # recalculate primitive vars (* wLC[3,:] is enthalpy now, not energy)
            uL = rhouL / rhoL
            uR = rhouR / rhoR
            vL = rhovR / rhoL
            vR = rhovR / rhoR
            pL = gmg*(rhohL - .5*(rhouL*rhouL + rhovL*rhovL)/rhoL)
            pR = gmg*(rhohR - .5*(rhouR*rhouR + rhovR*rhovR)/rhoR)

            # compute Roe average face quantities
            rhoLRt = sqrt(rhoL)
            rhoRRt = sqrt(rhoR)
            deNom = rhoLRt + rhoRRt
            rhoF = rhoLRt * rhoRRt;
            uF = (rhouL/rhoLRt + rhouR/rhoRRt) / deNom
            vF = (rhovL/rhoLRt + rhovR/rhoRRt) / deNom
            hF = (rhohL/rhoLRt + rhohR/rhoRRt) / deNom
            velMagSq = (uF**2 + vF**2)
            cFSq = (gm1)*(hF - .5*(velMagSq))
            cF = maximum(sqrt(cFSq),self.__cLim)

            normSF = sqrt(sum(sF*sF,0))
            qsL = sF[0,:]*uL + sF[1,:]*vL
            qsR = sF[0,:]*uR + sF[1,:]*vR
            qsF = sF[0,:]*uF + sF[1,:]*vF
            csF = cF*normSF
            qmF = qsF/csF
            #print cF.shape
            #print normSF.shape
            aa = (self.__dLim)*csF # should be 1 x nF

            rc = abs(qsF)
            ind = arange(rc.size)[rc < aa]
            rc[ind] = .5*(aa[ind] + ( (rc[ind]**2.)/aa[ind]) )
            #rp = sign(qmF)
            rp = ones(qmF.shape)
            rp[qmF < 0.] = -1.

            #print qsF.shape
            #print b.shape
            #print rp.shape
            #print rp.size

            ind = arange(qmF.size)[abs(qmF) < 1.]
            #print max(ind)
            #rp[ind] = maximum(2.*(abs(qmF[ind]))-1.,1e-6)*rp[ind] # flo75 Version
            ga = .5*(gma+1.)/gma
            bb = ga*abs(qsF[ind]) \
              - sqrt( (ga*qsF[ind])**2. + (csF[ind]**2. - qsF[ind]**2.)/gma )
            #print maximum(((abs(qsF[ind])+bb[ind])/(abs(qsF[ind])-bb[ind])),0.).shape
            rp[ind] = maximum(((abs(qsF[ind])+bb)/(abs(qsF[ind])-bb)),0.) \
                    * rp[ind]

            upLim = abs(self.__vis1)
            rc = upLim*(rc - rp*qsF)
            rp = upLim*rp        
            fc = rc + .5*rp*(qsL+qsR)
            fp = .5*rp*(qsR-qsL)
            fdp = rp*(pR-pL)

            fDisL = fc - fp
            fDisR = fc + fp

            fN = empty([nVar,msh.nE])
            fN[0,:] = fDisR*rhoR - fDisL*rhoL
            fN[1,:] = fDisR*rhouR - fDisL*rhouL + sF[0,:]*fdp
            fN[2,:] = fDisR*rhovR - fDisL*rhovL + sF[1,:]*fdp
            #fN[1,:] = fDisR*rhouR - fDisL*rhouL + sF[0,:]*fdp
            #fN[2,:] = fDisR*rhovR - fDisL*rhovL + sF[1,:]*fdp
            fN[3,:] = fDisR*rhohR - fDisL*rhohL

            #ind = arange(msh.nE)[e0 < nB]
            #fN[:,ind] = 0.
            #ind = arange(msh.nE)[e1 < nB]
            #fN[:,ind] = 0.

            fN = -fN
            #print 'fN Cusp'
            #for iV in range(nVar):
            #    for iF in range(msh.nE):
            #        print iV,':',iF,':', e0[iF], e1[iF], fN[iV,iF]

            #fNPT = zeros([nVar,msh.nP])
            #for iF in range(msh.nE):
            #    fNPT[:,e0[iF]] += fN[:,iF]
            #    fNPT[:,e1[iF]] -= fN[:,iF]
            #print 'fN Cusp @ point'
            #for iV in range(nVar):
            #    for iP in range(msh.nP):
            #        print iV,': point ',iP,':', fNPT[iV,iP]

            #sys.exit(0)
        return fN
    #---------------------------------------------------------------------------
    def getSelectedGrad(self,
                          wLC, wRC, pL, pR, # nF vectors
                          sF, msh,
                          w, p,
                          iF):

            nVar = wLC.shape[0]

            e0 = msh.e0[iF]
            e1 = msh.e1[iF]

            print list(hstack((e0,e1)))
            lP = list(set(list(hstack((e0,e1))))) # find unique points
            nP = len(lP)
            pMap0 = [lP.index(ii) for ii in e0]
            pMap1 = [lP.index(ii) for ii in e1]
            print lP

            ec0 = [msh.egClouds0[ii] for ii in lP]
            ec1 = [msh.egClouds1[ii] for ii in lP]

            print ec0, ec1

            ecS0 = hstack([ii for ii in ec0])
            ecS1 = hstack([ii for ii in ec1])
            print ecS0
            lF = list(set(list(hstack((ecS0,ecS1))))) # find unique points
            nF = len(lF)
            fMap0 = [lF.index(ii) for ii in ecS0]
            fMap1 = [lF.index(ii) for ii in ecS1]
            print lF
            print fMap0,fMap1
            
            sF = msh.sF[:,lF]

            sFe0 = sF/msh.mass[msh.e0[lF]]   # size of nD x nF
            sFe1 = sF/msh.mass[msh.e1[lF]]

            #WLC = w[:-1,e0].copy()
            #PL = w[4,e0].copy()
            #WRC = w[:-1,e1].copy()
            #PR = w[4,e1].copy()
            
            wLC[3,:] += pL  # wLC[3,:] is enthalpy
            wRC[3,:] += pR  # wRC[3,:] is enthalpy

            wF = .5*(wLC + wRC)

            wFNx0 = wF*tile(sFe0[0,:],[nVar,1])
            wFNy0 = wF*tile(sFe0[1,:],[nVar,1])

            wFNx1 = wF*tile(sFe1[0,:],[nVar,1])
            wFNy1 = wF*tile(sFe1[1,:],[nVar,1])

            wxC = empty([nVar,nP])
            wyC = empty([nVar,nP])
            for iV in range(nVar):
	        wxC[iV,:] = array([sum(jj,0)
                              for jj in
                                  zip([ sum(wFNx0[iV,ii]) for ii in ec0],
                                      [-sum(wFNx1[iV,ii]) for ii in ec1])
                              ])

	        wyC[iV,:] = array([sum(jj,0)
                              for jj in
                                  zip([ sum(wFNy0[iV,ii]) for ii in ec0],
                                      [-sum(wFNy1[iV,ii]) for ii in ec1])
                              ])

            # correct boundary values
            lB = [arange(nP)[ii < nB] for ii in e0]
            eB = e0[lB]
            pB = pMap0[lB]
            sB = msh.sB[:,eB]/tile(msh.mass[eB],[2,1])
            wxC[:,pB] += w[:-1,eB]*tile(sB[0,:],[nVar,1])
            wyC[:,pB] += w[:-1,eB]*tile(sB[1,:],[nVar,1])
            wxC[3,pB] += w[4,eB]*sB[0,:]
            wyC[3,pB] += w[4,eB]*sB[1,:]

            lB = [arange(nP)[ii < nB] for ii in e1]
            eB = e0[lB]
            pB = pMap1[lB]
            sB = msh.sB[:,eB]/tile(msh.mass[eB],[2,1])
            wxC[:,pB] += w[:-1,eB]*tile(sB[0,:],[nVar,1])
            wyC[:,pB] += w[:-1,eB]*tile(sB[1,:],[nVar,1])
            wxC[3,pB] += w[4,eB]*sB[0,:]
            wyC[3,pB] += w[4,eB]*sB[1,:]

            #wxC[:,:nB] += w[:-1,:nB]*tile(sB[0,:],[nVar,1])
            #wyC[:,:nB] += w[:-1,:nB]*tile(sB[1,:],[nVar,1])
            #wxC[3,:nB] += w[4,:nB]*sB[0,:]
            #wyC[3,:nB] += w[4,:nB]*sB[1,:]

            wx = wxC
            wy = wyC

            dx = msh.x[e1,:] - msh.x[e0,:]
            dwFL =   wx[:,pMap0]*(tile(dx[:,0],[nVar,1])) \
                   + wy[:,pMap0]*(tile(dx[:,1],[nVar,1]))
            dwFR =   wx[:,pMap1]*(tile(dx[:,0],[nVar,1])) \
                   + wy[:,pMap1]*(tile(dx[:,1],[nVar,1]))
            return dwFL,dwFR
##------------------------------------------------------------------------------
## end class InterfaceSolverCusp
##------------------------------------------------------------------------------
