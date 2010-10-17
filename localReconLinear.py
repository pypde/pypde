import sys
#import mesh
from numpy import *
#from numpy.linalg import *
#from pysparse import spmatrix

##------------------------------------------------------------------------------
class LocalReconLinear(object):
    #---------------------------------------------------------------------------
    # constructor
    #
    def __init__(self,eqn,msh,schmDict):

        print "     Linear reconstruction:"
        print "         Gradient:"
        for param,value in schmDict.get('gradient').asDict().iteritems():
            print "".join(["             ", param , ": ", str(value)])
        print "         Limiter:"
        for param,value in schmDict.get('limiter').asDict().iteritems():
            print "".join(["             ", param , ": ", str(value)])

                       
        self.nD = 2   # 2D

        gradDict = {}
        gradDict['gauss'] = self.buildGradGauss
        gradDict['leastSquares'] = self.buildGradLS
        gradDict['minNorm'] = self.buildGradMinNorm

        limiterDict = {}
        #limiterDict['alphaMean'] = slopeLimAlphaMean
	#limiterDict['barthJespersen'] = slopeLimBarthJespersen
        #limiterDict['minMod'] = slopeLimMinMod
        limiterDict['switchedAMean'] = slopeLimSwitchedAMean

        #gradArgDict = {}
        limArgDict = {}
        gradArgDict = {}

        gradFn = schmDict.get('gradient').get('type')
        limiter = schmDict.get('limiter').get('type')


        # limiter-specific argument dictionary
        switchedAMeanArgDict = {}
        switchedAMeanArgDict['eqn'] = eqn
        switchedAMeanArgDict['msh'] = msh

        limSpecArgDict = {}
        limSpecArgDict['switchedAMean'] = switchedAMeanArgDict



        if gradFn in gradDict:
            gradArgDict = schmDict.get('gradient').asDict()
            del gradArgDict['type']
            gradArgDict.update({'msh':msh})
            self.__grad = gradDict.get(gradFn)(msh=msh)
            if len(gradArgDict) >0:
                self.slopeLimiter = gradDict.get(gradFn)(**gradArgDict)

            else:
                self.slopeLimiter = gradDict.get(gradFn)()
        else:
            sys.exit(' *** ERROR: Unsupported gradient scheme')


        if limiter in limiterDict:
            limArgDict = schmDict.get('limiter').asDict()
            limArgDict.update(limSpecArgDict.get(limiter))
            del limArgDict['type']
            if len(limArgDict) >0:
                self.slopeLimiter = limiterDict.get(limiter)(**limArgDict)
    
            else:
                self.slopeLimiter = limiterDict.get(limiter)()
        else:
            sys.exit(' *** ERROR: Unsupported gradient limiter')

    #---------------------------------------------------------------------------
    def buildGradGauss(self,msh):
        return self.__gradGauss

    def __gradGauss(self,w,msh):
       
        nVar = w.shape[0]
        gradW = zeros([self.nD,nVar,msh.nP])   # local nodal solution gradient
        sF = msh.sF
        nE = msh.nE
        nB = msh.nB
        ec0 = msh.egClouds0
        ec1 = msh.egClouds1
        e0 = msh.e0
        e1 = msh.e1

        sFe0 = sF/msh.mass[e0]
        sFe1 = sF/msh.mass[e1]

        # use mesh to obtain gradient
        for e in range(nE):
            gradW[0,:,e0[e]] += w[:,e1[e]]*sFe0[0,e]
            gradW[0,:,e1[e]] -= w[:,e0[e]]*sFe1[0,e]
        for e in range(nE):
            gradW[1,:,e0[e]] += w[:,e1[e]]*sFe0[1,e]
            gradW[1,:,e1[e]] -= w[:,e0[e]]*sFe1[1,e]

        # correct boundary values (need to get original w somehow
        sB = msh.sB/tile(msh.mass[:nB],[2,1])
        gradW[0,:,:nB] += w[:,:nB]*tile(sB[0,:],[nVar,1])
        gradW[1,:,:nB] += w[:,:nB]*tile(sB[1,:],[nVar,1])

        return .5*gradW # (b/c meshless aij is sF/2)
    #---------------------------------------------------------------------------

    def buildGradLS(self,msh,weighting='none'):
    #---------------------------------------------------------------------------
    # build local linear derivative operators
    #---------------------------------------------------------------------------
    # least squares (linear, assuming more than 2 neighbors)

        self.bij = []
        for i in range(msh.nP):
            #nbrInd = arange(msh.nbrS[i],msh.nbrS[i+1])
            dx = msh.x[hstack([msh.clouds0[i],msh.clouds1[i]]),:] \
               - tile(msh.x[i,:],[msh.nNbr[i],1])
            #print dx
            #dx = msh.x[msh.clouds[nbrInd,1],:] - msh.x[msh.clouds[nbrInd,0],:]

            if weighting == 'none':
                weights = ones([msh.nNbr[i],1],float64)
            elif weighting == 'inverseDist':
                weights = 1. / sqrt(sum(dx**2,1))
            elif weighting == 'inverseDistPower2':
                weights = 1. / sum(dx**2,1)
            else:
                sys.exit(' *** ERROR: Unsupported weighting scheme for \
least squares gradient. ***')

            #wx = sum(weights * dx[:,0],1)
            #wy = sum(weights * dx[:,1],1)
            wxy = sum(weights * dx[:,0]* dx[:,1])
            wx2 = sum(weights * dx[:,0]**2)
            wy2 = sum(weights * dx[:,1]**2)

            denom = wx2 * wy2 - wxy**2
            bijx =  (weights*dx[:,0]*wy2 - weights*dx[:,1]*wxy) / denom
            bijy =  (weights*dx[:,1]*wx2 - weights*dx[:,0]*wxy) / denom
            #print bijx
            #print bijy
            #print vstack([bijx,bijy]).transpose()
            # [b_{ij}^1, b_{ij}^2] for local derivatives, (nNbr x 2)
            self.bij.append(vstack([bijx,bijy]).transpose())

        return self.__gradLS
    #---------------------------------------------------------------------------
    def __gradLS(self,w,msh):

        nVar = w.shape[0]
        gradW = zeros([self.nD,nVar,msh.nP])   # local nodal solution gradient
        # cell to edge extrapolation using least squares
        for id in range(self.nD):
            for iS in range(msh.nP):
                #nbrInd = msh.clouds[arange(msh.nbrS[iS],msh.nbrS[iS+1]),1]
                nbrInd = hstack([msh.clouds0[iS],msh.clouds1[iS]])
                gradW[id,:,iS] = \
                    sum(  
                          tile(array(self.bij[iS])[:,id],[nVar,1]) \
                        * (  w[:,nbrInd] \
                           - tile(w[:,iS],[msh.nNbr[iS],1]).transpose()),
                        1)

        return gradW
    #---------------------------------------------------------------------------
    def buildGradMinNorm(self,msh):
        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # build local linear derivative operators
        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # min norm
        self.bij = []
        for i in range(msh.nP):
            nbrInd = arange(msh.nbrS[i],msh.nbrS[i+1])
            dx = msh.x[msh.clouds[nbrInd,1],:] - msh.x[msh.clouds[nbrInd,0],:]
            C = hstack([ones([msh.nNbr[i],1]), dx])
            Q, R = qr(C)
            # a = dot(Q, solve(R.transpose(), [[0,0],[1,0],[0,1]]))
            # gradOp = egClouds[nbrInd,1]*(a * dx).sum(1)

            # [b_{ij}^1, b_{ij}^2] for local derivatives, (nNbr x 2)
            self.bij.append(dot(Q, solve(R.transpose(), [[0,0],[1,0],[0,1]])))  
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    def sp2Fp (self, w, msh): #,nSym='',eSym=''):

        nVar = w.shape[0]

        dxe0 = msh.xE - msh.x[msh.e0,:] # edge vectors
        dxe1 = msh.xE - msh.x[msh.e1,:] # edge vectors

        gradW = self.__grad(w,msh)

        dwL = sum(
                   gradW[:,:,msh.e0] \
                 * tile((dxe0.transpose())[:,newaxis,:],[1,nVar,1]),
                 0)

        #print dwL.shape

        dwR = sum(
                    gradW[:,:,msh.e1] \
                  * tile((dxe1.transpose())[:,newaxis,:],[1,nVar,1]),
                  0)

        # CAN & SHOULD change limiter function to return the LIMITED DW BY EDGE,
        # not just the limiter!
        #for iV in range(nVar):
        #    phi[iV,:] = self.slopeLimiter.limitSlope(w[iV,:],dwL[iV,:],dwR[iV,:],msh)
        #phi = tile(self.slopeLimiter.limitSlope(w[0,:],dwL[0,:],dwR[0,:],msh),[5,1])

        #wL = w[:,msh.e0] + phi[:,msh.e0]*dwL
        #wR = w[:,msh.e1] + phi[:,msh.e1]*dwR

        dwLLim, dwRLim = self.slopeLimiter.limitSlope(w,dwL,dwR,msh)

        wL = w[:,msh.e0] + dwLLim
        wR = w[:,msh.e1] + dwRLim

        return wL, wR

    #---------------------------------------------------------------------------
##------------------------------------------------------------------------------
## end class LocalReconLinear
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
class slopeLimSwitchedAMean(object):
    def __init__ (self, eqn=None, msh=None, q=3., q2=5., vis=1., wLim=1e-6):
        self.__q = q # q=1:minmod; q=2:van Leer for q2 = 1.
        self.__q2 = q2
        self.__wLim = wLim
        #bLim = float(vis)/(float(msh.nP)**(.75))
        #rho0 = eqn.refState().rho()
        #c0 = eqn.refState().c()
        #h0 = eqn.refState().h()
        #p0 = eqn.refState().p()
        #self.__wLimVec = bLim*array([rho0,rho0*c0,rho0*c0,rho0*h0-p0])
        #self.__wLim = bLim*p0
        self.__wLim = wLim

    def limitSlope(self, w, dwL, dwR, msh):

        nVar = dwL.shape[0]
        nE = dwL.shape[1]
	nWSize = dwL.size

        
        #if nVar == 1:
        #    # scalar limit coefficient
        #    wLimA = tile(self.__wLim,[nVar,nE])
        #elif nVar == 4:
        #    # vector limit coefficients
        #    wLimA = tile(self.__wLimVec,[nE,1]).transpose()
        #else:
        #    sys.error(' *** Error: unexpected number of variables for limiter!')
        wLimA = self.__wLim

        aShape = [nVar,nE]
        a = reshape(maximum(abs(dwL) + abs(dwR),wLimA),nWSize)
        b = reshape(abs(dwR + dwL),nWSize) # note signage!
        r = ones(nWSize)

        #ind = arange(r.size)[a > reshape(wLimA,nWSize)]
        ind = arange(r.size)[a > wLimA]
        r[ind] = 1. - (b[ind]/a[ind])**self.__q # q->0, dissipation->1st order.a
        #dwLim = .5*reshape(r,aShape)*(dwL - dwR) # limited average of dw's
        dwLim = .5*reshape(r**(self.__q2),aShape)*(dwL - dwR) # limited average of dw's

        return dwLim, -dwLim
    
###------------------------------------------------------------------------------
#class slopeLimAlphaMean(object):
#    def __init__ (self,alpha):
#        self.alpha = alpha
#     
#    def limitSlope(self, dwL, dwR):
#
#        s = .5*(sign(dwL) + sign(dwR))
#        phi = s * amin( [abs(dwL + dwR)/2, 
#                       abs(dwL)*self.alpha, 
#                       abs(dwR)*self.alpha],0)
#        return phi
#    
###------------------------------------------------------------------------------
#class slopeLimMinMod(object):
#     
#    def limitSlope(self, dwL, dwR):
#
#        s = .5*(sign(dwL) + sign(dwR))
#        phi = s * amin([abs(dwL),abs(dwR)],0)
#        return phi
#
###------------------------------------------------------------------------------
class slopeLimBarthJespersen(object):
     
    def limitSlope(self, w, dwL, dwR, msh):

        #tol = finfo(double).eps
        tol = 1.e-6
        phi = zeros(msh.nP)

        for i in range(msh.nP):
            #phiE = ones(msh.nNbr[i])
            phiE0 = ones(msh.nNbrP[i])
            phiE1 = ones(msh.nNbrM[i])
            #nbrInd = arange(msh.nbrS[i], msh.nbrS[i+1])
            #print w[msh.clouds[nbrInd,1]].shape
            wLocal = hstack([w[msh.clouds0[i]],w[msh.clouds1[i]],w[i]])
            wMax  = max(wLocal)
            wMin  = min(wLocal)

            ## local indices of neighbors such that i is e0 of the edge
            #nbrInd0 = arange(msh.nNbr[i])[msh.egClouds[nbrInd,1]== 1] 
            ## local indices of neighbors such that i is e1 of the edge
            #nbrInd1 = arange(msh.nNbr[i])[msh.egClouds[nbrInd,1]==-1]

            ## global indices of edges such that i is e0 of the edge
            #egInd0 = msh.egClouds[ nbrInd0 + msh.nbrS[i] , 0]
            ## global indices of edges such that i is e1 of the edge
            #egInd1 = msh.egClouds[ nbrInd1 + msh.nbrS[i] , 0]

            #indP0, indP1, indM0, indM1 = empty(0),empty(0),empty(0),empty(0)
            if (abs(wMax-w[i]) > tol):
                indP0 = arange(msh.nNbr[i],dtype=int) \
                            [dwL[msh.egClouds0[i]] >  tol]
                indP1 = arange(msh.nNbr[i],dtype=int) \
                            [dwR[msh.egClouds1[i]] >  tol]
            else:
                indP0, indP1 = empty(0,dtype=int),empty(0,dtype=int)
              
            if (abs(wMin-w[i]) > tol):
                indM0 = arange(msh.nNbr[i],dtype=int) \
                            [dwL[msh.egClouds0[i]] < -tol]
                indM1 = arange(msh.nNbr[i],dtype=int) \
                            [dwR[msh.egClouds1[i]] < -tol]
            else:
                indM0, indM1 = empty(0,dtype=int),empty(0,dtype=int)

            if indP0.shape[0] > 0:
            #if (abs(wMax-w[i])>tol) and indP0.shape[0] > 0:
               #phiE[nbrInd0[indP0]] = (wMax-w[i])/dwL[msh.egClouds0[i][indP0]]
               phiE0[indP0] = (wMax-w[i])/dwL[msh.egClouds0[i][indP0]]
            if indM0.shape[0] > 0:
            #if (abs(wMin-w[i])>tol) and indM0.shape[0] > 0:
               #phiE[nbrInd0[indM0]] = (wMin-w[i])/dwL[msh.egClouds0[i][indM0]]
               phiE0[indM0] = (wMin-w[i])/dwL[msh.egClouds0[i][indM0]]
            if indP1.shape[0] > 0:
            #if (abs(wMax-w[i])>tol) and indP1.shape[0] > 0:
               #phiE[nbrInd1[indP1]] = (wMax-w[i])/dwR[msh.egClouds1[i][indP1]]
               phiE1[indP1] = (wMax-w[i])/dwR[msh.egClouds1[i][indP1]]
            if indM1.shape[0] > 0:
            #if (abs(wMin-w[i])>tol) and indM1.shape[0] > 0:
               #phiE[nbrInd1[indM1]] = (wMin-w[i])/dwR[msh.egClouds1[i][indM1]]
               phiE1[indM1] = (wMin-w[i])/dwR[msh.egClouds1[i][indM1]]

            #if  (i == 9 or i == 10):
            #    print 'i:', i
            #    print 'wMax:', wMax
            #    print 'wMin:', wMin
            #    print indP0, indM0, indP1, indM1
            #    print indP0.shape, indM0.shape, indP1.shape, indM1.shape
            #    print 'dw''s:', \
            #          dwL[msh.egClouds0[i][indP0]], \
            #          dwL[msh.egClouds0[i][indM0]], \
            #          dwR[msh.egClouds1[i][indP1]], \
            #          dwR[msh.egClouds1[i][indM1]]
            #    print phiE0
            #    print phiE1

            #print phiE
            #print hstack([1., phiE])
            #phi[i] = min(hstack([1., phiE]))

            phi[i] = min(hstack([1., phiE0, phiE1]))

        return phi
