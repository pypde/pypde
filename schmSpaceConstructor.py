import sys

##------------------------------------------------------------------------------
def constructSchmSpace(eqn,msh,spaceDict):

    reconDict = dict(spaceDict.get('reconstruction'))

    if spaceDict.get('family') == 'finiteVolume':

        from schmSpaceFVFluxRecon import SchmSpaceFVFluxRecon
        from schmSpaceFVSolRecon import SchmSpaceFVSolRecon

        schmFVConstructorDict = {}
        schmFVConstructorDict['flux'] = SchmSpaceFVFluxRecon
        schmFVConstructorDict['solution'] = SchmSpaceFVSolRecon

        # sort out whether to reconstruct solution or flux
        reconQty = reconDict.get('quantity')
        reconArgs = reconDict.copy()
        del reconArgs['quantity']

        schmSpace = schmFVConstructorDict[reconQty](eqn,msh,reconArgs)

        #if (reconQty == 'solution'):
        #    sys.exit(' *** ERROR: Solution reconstruction \
        #                           currently unsupported. ***')
        #elif (reconQty == 'flux'):
        #    schmSpace = SchmSpaceFVFluxRecon(msh,reconArgs)
    else:
         sys.exit(' *** ERROR: Unsupported spatial scheme! ***')

    return schmSpace
##------------------------------------------------------------------------------
