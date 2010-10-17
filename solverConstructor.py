#from explicitRK import *

##------------------------------------------------------------------------------
def constructSolver(solverDict):

    solverModDict = {}
    solverModDict['explicitRK'] = 'solverExplicitRK'
    solverModDict['jamesonRK'] = 'solverJamesonRK'

    solverConstructorDict = {}
    solverConstructorDict['explicitRK'] = 'SolverExplicitRK'
    solverConstructorDict['jamesonRK'] = 'SolverJamesonRK'

    solverArgs = solverDict.copy()
    del solverArgs['family']

    #timeStepMode = solverDict.get('timeStepMode')
    #if timeStepMode == 'local':
    #    tModeLocal = 1
    #elif timeStepMode == 'global':
    #    tModeLocal = 0
    #else:
    #    sys.exit(' *** ERROR: Unsupported time step mode. ***')

    #del solverArgs['timeStepMode']
    #solverArgs.update({'tModeLocal':tModeLocal})


    if solverDict['family'] in solverConstructorDict:
        ts = __import__(solverModDict.get(solverDict.get('family')))
        solver = getattr(ts,solverConstructorDict[solverDict['family']])(**solverArgs)
    else:
        sys.exit(' *** ERROR: Unsupported temporal scheme / solver! ***')

    #solver = solverConstructorDict[solverDict['family']]\
    #                (**solverArgs)

    return solver
##------------------------------------------------------------------------------
