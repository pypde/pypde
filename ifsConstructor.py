from copy import deepcopy
from ifsCentral2 import InterfaceSolverCentral2
from ifsCusp import InterfaceSolverCusp
from ifsRoe import InterfaceSolverRoe

##------------------------------------------------------------------------------
def constructInterfaceSolver(eqn,msh,ifsDict):

    ###-------------------------------------------------------------------------
    ### set up lookup tables for class constructors
    ###-------------------------------------------------------------------------
    ifsConstructorDict = {}
    ifsConstructorDict['central2'] = InterfaceSolverCentral2
    ifsConstructorDict['cusp'] = InterfaceSolverCusp
    ifsConstructorDict['roe'] = InterfaceSolverRoe

    # list common arguments to pass
    ifsConstructorSharedArgDict = {}
    ifsConstructorSharedArgDict['msh'] = msh

    # list ifs specific arguments (can add ifs-specific ones here)
    ifsConstructorArgDict = {}

    ifsCuspArgDict = {}
    ifsCuspArgDict['refState'] = eqn.refState()

    ifsRoeArgDict = {}
    ifsRoeArgDict['gma'] = eqn.gma()
    
    ifsConstructorArgDict['central2'] = {}
    ifsConstructorArgDict['cusp'] = ifsCuspArgDict
    ifsConstructorArgDict['roe'] = ifsRoeArgDict

    ifss = []
    for ifsType in (ifsDict.keys()):

        # get parameters for interface solver based on input file
        thisIfsArgs = dict(ifsDict.get(ifsType))
        # get extra arguments specific to the current interface solver 
        thisIfsArgs.update(ifsConstructorArgDict.get(ifsType))
        # get extra arguments applicable for all interface solvers
        thisIfsArgs.update(ifsConstructorSharedArgDict)

        #print thisIfsArgs
        if len(thisIfsArgs) > 0:
            #print ifsConstructorDict[ifsType]
            ifss.append(ifsConstructorDict[ifsType](**thisIfsArgs))
        else:
            ifss.append(ifsConstructorDict[ifsType]())

    #print ifss
    return ifss
##------------------------------------------------------------------------------
