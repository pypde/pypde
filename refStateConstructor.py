from refStateFreeStream import *

##------------------------------------------------------------------------------
def constructRefState(eqn,refStateDict):

    refStateConstructorDict = {}
    refStateConstructorDict['freeStream'] = RefStateFreeStream

    refStateArgs = refStateDict.copy()
    del refStateArgs['type']

    refState = refStateConstructorDict[refStateDict['type']]\
                    (eqn,**refStateArgs)

    return refState
##------------------------------------------------------------------------------
