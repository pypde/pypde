from localReconConstant import *
from localReconLinear import *

##------------------------------------------------------------------------------
def constructSchmLocalRecon(eqn,msh,reconSchmDict):

    reconSchmConstructorDict = {}
    reconSchmConstructorDict['constant'] = LocalReconConstant
    reconSchmConstructorDict['linear'] = LocalReconLinear

    reconSchmArgs = reconSchmDict.copy()
    del reconSchmArgs['accuracy']

    reconSchm = reconSchmConstructorDict[reconSchmDict['accuracy']]\
                    (eqn,msh,reconSchmArgs)

    return reconSchm
##------------------------------------------------------------------------------
