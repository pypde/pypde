import sys
import operator
from pyparsing import *
from numpy import *


class odict(dict):
    def __init__(self, *argt, **argd):
        dict.__init__(self, *argt, **argd)
        self.__dict__ = self


##------------------------------------------------------------------------------
class Mesh(object):
    def __init__(self,filenameMesh):

        ##......................................................................
        ## Parse mesh file
        ##......................................................................
        
        print ' ... Parsing mesh file'
        sys.stdout.flush()
        fM = open(filenameMesh,'r')
        expressions = fM.read();

        EQ,LBRACE,RBRACE,SEMI = map(Suppress,'={};')
        LPAREN,RPAREN = map(Suppress,'()')

        intEntry = Combine(Optional('-') + Word(nums)) \
                   .setParseAction(lambda t : int(t[0]))
        realEntry = ( Combine(Optional('-') + Word(nums) \
                            + '.' + Optional(Word(nums))) \
                    | Combine(Optional('-') + Optional(Word(nums)) \
                            + '.' + Word(nums)) )
        numEntry = realEntry | intEntry
        realOrSciEntry = Combine(  numEntry \
                                 + Optional(CaselessLiteral('E') \
                                            + Word('+-'+nums,nums))) \
                         .setParseAction(lambda t : float(t[0]))

       # define tokens, expressions and entries
        objEntry = Forward()
        vecEntry = Forward()

        keyToken = Word(alphas+'_', alphanums+'_')
        entryToken = ( keyToken | intEntry | realOrSciEntry |
                       quotedString.copy().setParseAction(removeQuotes))

        # define lists
        numList = Group(OneOrMore(realOrSciEntry))
        intList = Group(OneOrMore(intEntry))
        #keyTokenList = Group(OneOrMore(keyToken))
        tupleList = (  Group(OneOrMore(LPAREN + intList + RPAREN)) \
                     | Group(OneOrMore(LPAREN + numList + RPAREN)))

        objExpr = (Group(keyToken + objEntry)|Group(keyToken + LBRACE + RBRACE))
        #vecExpr =  Group(intEntry + vecEntry)
        vecExpr =  (intEntry + vecEntry)
        expr    =  Group(keyToken + entryToken + SEMI)

        genExpr = (expr | objExpr)
        mixedExpr = (vecExpr | expr | objExpr)

        vecEntry << (  (LPAREN + intList + RPAREN) \
                     | (LPAREN + numList + RPAREN) \
                     | (LPAREN + tupleList + RPAREN))
        objEntry << (  (  LBRACE + Dict(OneOrMore(genExpr)) + RBRACE) \
                     | (  LBRACE + OneOrMore(vecExpr) + RBRACE))

        mshDict = Dict(OneOrMore(mixedExpr)) \
                       .ignore(cStyleComment).ignore('//' + restOfLine) \
                       .parseString(expressions,parseAll=True)
        
        print '     Done'
        sys.stdout.flush()
        #print mshDict.dump()

        print ' ... Building domain'
        sys.stdout.flush()
        nP = mshDict['nodeCoordinates'][0]
        assert (nP == mshDict['volumes'][0])

        nB = mshDict['boundaryFaceAreas'][0]

        nE = mshDict['edgeConnectivity'][0]
        assert (nE == mshDict['interiorFaceCenters'][0])
        assert (nE == mshDict['interiorFaceAreas'][0])

        assert (nP == len(mshDict['nodeCoordinates'][1]))
        assert (nP == len(mshDict['volumes'][1]))
        assert (nB == len(mshDict['boundaryFaceAreas'][1]))
        assert (nE == len(mshDict['interiorFaceCenters'][1]))
        assert (nE == len(mshDict['interiorFaceAreas'][1]))
        assert (nE == len(mshDict['edgeConnectivity'][1]))

        x = array(mshDict['nodeCoordinates'][1].asList())
        sB = array(mshDict['boundaryFaceAreas'][1].asList()).transpose()
        mass = array(mshDict['volumes'][1].asList())

        xE = array(mshDict['interiorFaceCenters'][1].asList())
        sF = array(mshDict['interiorFaceAreas'][1].asList()).transpose()
        #print (mshDict['interiorFaceAreas'][1].asList())
        #print sF.dtype
        #print mass.dtype

        #print 'sF'
        #print sF.transpose()
        #sys.exit(0)
        #print 'sB'
        #print sB.transpose()
        #print 'xE'
        #print xE
        #print 'x'
        #print x

        e0 = array(mshDict['edgeConnectivity'][1].asList())
        e1 = e0[:,1]-1
        e0 = e0[:,0]-1
        #e0 = (mshDict['edgeConnectivity'][1].asList())
        #e1 = [i[1]-1 for i in e0]
        #e0 = [i[0]-1 for i in e0]

        # check and make sure mesh is oriented the right way...?
        #flipIndices = arange(nE)[sum(sF.transpose()*(x[e1,:] - x[e0,:]),1) < 0]
        #e0[flipIndices],e1[flipIndices] = e1[flipIndices],e0[flipIndices]
        #print ' *** flipped ', flipIndices.shape[0], ' edges.'

        #ttt = zip(arange(nE), e0, e1)
        #for iii in range(nE):
        #    print ttt[iii]
        #print e0.shape ,e1.shape

        bcDict_t = mshDict['boundaryPatches']
        assert( nB == sum(bcDict_t[iB]['nFaces'] for iB in bcDict_t.keys()))
        #assert( nB == sum(patchList[iB][1] for iB in range(len(patchList))))
        patchList = []
        for i in bcDict_t.keys():
            patchList.append((i,bcDict_t[i]['nFaces'],bcDict_t[i]['startFace']-1))  # convert from 1-based to 0-based
        
        patchList.sort(key=operator.itemgetter(2))
        #print patchList

        # number of boundary patches
        nPatches = len(patchList)
        # bounding box
        xBB = vstack([[min(x[:,0]),min(x[:,1])],[max(x[:,0]),max(x[:,1])]])

        ##......................................................................
        ## print point statistics
        ##......................................................................
        print '     done reading points and edges.'
        print '     number of total points:         ', nP
        #print '     number of interior points:      ', nPI
        print '     number of boundary points:      ', nB
        print '     number of patches:              ', nPatches
        print '     number of points per patch:'
        for i in range(len(patchList)):
            print '        ', patchList[i][0], ': ', \
                              patchList[i][1]
        print '     number of total edges:   ', nE
        print '     bounding box:            ', xBB[0,:] ,xBB[1,:]
        sys.stdout.flush()

        # read in the left and right neighbors of the edge
        #clouds = []   # holds clouds for interior points in stacked from
        
        print ' ... Building info for neighbor connectivity'

        nNbrP = zeros(nP)
        nNbrM = zeros(nP)
        nNbrP[:max(e0)+1] = bincount(e0)
        nNbrM[:max(e1)+1] = bincount(e1)
        nNbr = nNbrP + nNbrM
        #print nNbrP, nNbrM
        #print nNbrP.shape, nNbrM.shape
        #print nNbr

        ##for i in range(nP):
        ##    print sorted(hstack([e1[e0==i],e0[e1==i]]))
        ##    print [e0[e1[j]==i for j in range(nP)]]
        ##clouds = hstack([vstack(( i*ones(nNbr[i],dtype=int),
        ##                          sorted(hstack([e1[e0==i],e0[e1==i]])) ) 
        ##                       )
        ##                 for i in range(nP)])
        #clouds = hstack([vstack(( i*ones(nNbr[i],dtype=int),
        #                          hstack([e1[e0==i],e0[e1==i]]) ) 
        #                       )
        #                 for i in range(nP)])
        ##print clouds

        #egClouds = hstack([ hstack([arange(nE)[e0==i],
        #                            arange(nE)[e1==i]])
        #                   for i in range(nP)])
        ##egClouds = hstack([ vstack(( hstack([arange(nE)[e0==i],
        ##                                     arange(nE)[e1==i]]),
        ##                             hstack([ ones(nNbrP[i],dtype=int),
        ##                                     -ones(nNbrM[i],dtype=int)])))
        ##                  for i in range(nP)])
        ##print egClouds

        clouds0 = [e1[e0==i] for i in range(nP)]
        egClouds0 = [arange(nE)[e0==i] for i in range(nP)]
        clouds1 = [e0[e1==i] for i in range(nP)]
        egClouds1 = [arange(nE)[e1==i] for i in range(nP)]
        #print 'clouds0'
        #print clouds0
        #print 'clouds1'
        #print clouds1
        #print 'egClouds0'
        #print egClouds0
        #print 'egClouds1'
        #print egClouds1


        #for i in range(nE):
        #    line = fM.readline().split()
        #    e0.append(int(line[0])-1)   # convert from 1-based to 0-based
        #    e1.append(int(line[1])-1)
        #    xE.append([float(line[2]), float(line[3])])
        #    sF.append([float(line[4]), float(line[5])])
        #    if (sum(sF[-1]*(x[e1[-1],:] - x[e0[-1],:])) < 0):
        #       print('flipped edge %d\n',i)
        #       e0[-1],e1[-1] = e1[-1],e0[-1] 
        #    clouds.append([e0[-1],e1[-1],i,1])
        #    clouds.append([e1[-1],e0[-1],i,-1])
        #
        #sorted(clouds,key=operator.itemgetter(1))
        #clouds.sort()
        #clouds = array(clouds)
        #egClouds = clouds[:,2:]
        #clouds = clouds[:,:2]
        #sF = array(sF).transpose()
        #xE = array(xE)

        #print e0
        #print e1
        #print clouds
        #print egClouds
        #
        ##get number of neighbors for each cloud
        #nNbr = []
        #for ll in range(nP):
        #    nNbr.append(e0.count(ll)+e1.count(ll))
        #
        ##build indexing and neighbor arrays
        #nbrS = reshape(hstack((0, cumsum(nNbr))),(nP+1,1))
        #nbr = []
        #for ll in range(nP):
        #    nbr.append(clouds[nbrS[ll]:nbrS[ll+1],1])
        #
        #e0, e1 = array(e0), array(e1)
        #
        ##print egClouds
        ##print clouds
        ##print e0
        ##print e1
        #print nbr
        #print nbrS
         
        print '     Done'
        sys.stdout.flush()
        fM.close()

        self.nP = nP
        self.nE = nE
        self.nB = nB
        self.nPatches = nPatches
        self.BB = xBB  # bounding box
        self.patchList = patchList

        self.xE = xE
        self.x = x
        self.mass = mass
        self.sF = sF
        self.sB = sB
        self.e0 = e0
        self.e1 = e1
        self.clouds0 = clouds0
        self.clouds1 = clouds1
        self.egClouds0 = egClouds0
        self.egClouds1 = egClouds1

        #self.clouds = clouds
        #self.egClouds = egClouds
        self.nNbr = nNbr
        self.nNbrP = nNbrP
        self.nNbrM = nNbrM
        #self.nbrS = nbrS
    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    def nP(self):
        return self.nP
    ##--------------------------------------------------------------------------
    ##--------------------------------------------------------------------------
    def numPatches(self):
        return self.nPatches
    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    def patchNum(self,patchName):
        nameList = [bcItem[0] for bcItem in self.patchList]
        return (nameList).index(patchName)
    ##--------------------------------------------------------------------------
    
    ##--------------------------------------------------------------------------
    def patchName(self,patchNum):
        return self.patchList[patchNum][0]
    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    def numFacesOnPatch(self,patchNum):
        return self.patchList[patchNum][1]
    ##--------------------------------------------------------------------------
    
    ##--------------------------------------------------------------------------
    def startFaceOfPatch(self,patchNum):
        return self.patchList[patchNum][2]
    ##--------------------------------------------------------------------------
    
    ##--------------------------------------------------------------------------
    def endFaceOfPatch(self,patchNum):
        if (patchNum < self.numPatches()-1):
             return self.patchList[patchNum+1][2]
        else:
             return self.nB
    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    def facesOnPatch(self,patchNum):
        return arange(self.startFaceOfPatch(patchNum),self.endFaceOfPatch(patchNum))
    ##--------------------------------------------------------------------------
    


