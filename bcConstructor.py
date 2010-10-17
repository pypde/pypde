#from copy import deepcopy
#from bcClFarfieldCirc import *
#from bcFarfieldCirc import *
#from bcFarfieldConst import *
#from bcFarfieldDoublet import *
#from bcVelTangentWall import *
import operator
import sys

##------------------------------------------------------------------------------
#def constructBC(eqn,msh,filenameBC):
def constructBC(eqn,refState,msh,bcDict):

    ##--------------------------------------------------------------------------
    ## set up lookup tables for class constructors
    ##--------------------------------------------------------------------------
    bcModuleDict = {}

    # module names
    eulerBCModuleDict = {}
    eulerBCModuleDict['fixedLiftCirculationCorrection'] = 'bcClFarfieldCirc'
    eulerBCModuleDict['circulationCorrection'] = 'bcFarfieldCirc'
    eulerBCModuleDict['constant'] = 'bcFarfieldConst'
    eulerBCModuleDict['doubletCorrection'] = 'bcFarfieldDoublet'
    eulerBCModuleDict['wall'] = 'bcVelTangentWall'

    bcModuleDict['Euler'] = eulerBCModuleDict

    bcConstructorDict = {}
    # euler equations
    eulerBCConstructorDict = {}
    eulerBCConstructorDict['fixedLiftCirculationCorrection'] = 'BCClFarfieldCirc'
    eulerBCConstructorDict['circulationCorrection'] = 'BCFarfieldCirc'
    eulerBCConstructorDict['constant'] = 'BCFarfieldConst'
    eulerBCConstructorDict['doubletCorrection'] = 'BCFarfieldDoublet'
    eulerBCConstructorDict['wall'] = 'BCVelTangentWall'

    bcConstructorDict['Euler'] = eulerBCConstructorDict

    # list common arguments to pass
    bcConstructorSharedArgDict = {}
    bcConstructorSharedArgDict['eqn'] = eqn
    bcConstructorSharedArgDict['msh'] = msh
    

    # list bc specific arguments (can add equation specific ones here)
    bcConstructorArgDict = {}
    eulerBCFarfieldArgDict = {}
    eulerBCFarfieldArgDict['refState'] = refState

    eulerVelTngtWallArgDict = {}

    bcConstructorArgDict['fixedLiftCirculationCorrection'] = eulerBCFarfieldArgDict
    bcConstructorArgDict['circulationCorrection'] = eulerBCFarfieldArgDict
    bcConstructorArgDict['constant'] = eulerBCFarfieldArgDict
    bcConstructorArgDict['doubletCorrection'] = eulerBCFarfieldArgDict
    bcConstructorArgDict['wall'] = eulerVelTngtWallArgDict

    # list bc specific arguments (can add equation specific ones here)
    for iDict in bcConstructorArgDict.values():
        iDict.update(bcConstructorSharedArgDict)
 
    # make a list of the boundary condition type (in order of patch number)
    bcTypeList = []
    for i in range(len(msh.patchList)):
        bcTypeList.append(bcDict[msh.patchList[i][0]]['type'])

    bcs = []
    for iB in range(len(bcTypeList)):
    #for iB in range(len(bcTypeList)-1,0,-1):

        if bcTypeList[iB] in bcConstructorDict.get(eqn.type):
            bcMod = __import__(bcModuleDict.get(eqn.type) \
                               .get(bcTypeList[iB]))

            bcName = msh.patchList[iB][0]

            # form arguments to pass to constructor for the i-th boundary
            thisBCArgs = {}

            if 'interfaceSolvers' in bcDict[bcName].keys():
                ifsDict = bcDict[bcName].get('interfaceSolvers').asDict()
                thisBCArgs.update(dict({'ifsDict':ifsDict}))
               
            for iEnt in (bcDict[bcName].keys()):
                if iEnt != 'interfaceSolvers':
                    thisBCArgs.update( dict({iEnt:bcDict[bcName].get(iEnt)})) # actual copy, not aliased
       
            del thisBCArgs['type']
            thisBCArgs['patchNum'] = iB
            thisBCArgs.update(bcConstructorArgDict[bcTypeList[iB]])
            thisBCArgs.update(bcConstructorSharedArgDict)

            #print thisBCArgs
            #bcs.append(bcConstructorDict[eqn.type][bcTypeList[iB]]\
            #                        (**thisBCArgs))
            bcs.append(getattr(bcMod,bcConstructorDict.get(eqn.type)\
                               .get(bcTypeList[iB]))\
                       (**thisBCArgs))
        else:
            sys.exit(' *** ERROR: Unsupported boundary condition! ***')


    # sort bcs according to priority
    print 'before sorting:'
    print bcs
    bcs.sort(key=operator.attrgetter('priority'))

    print 'after sorting:'
    print bcs
    return bcs
##------------------------------------------------------------------------------
