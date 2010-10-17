##------------------------------------------------------------------------------
def constructFigures(eqn,msh,figureDict,prefDir):

    ###--------------------------------------------------------------------------
    ### set up lookup tables for class constructors
    ###--------------------------------------------------------------------------
    figureConstructorDict = {}
    figureConstructorDict['mlab'] = 'FigureMlab'

    figureModuleDict = {}
    figureModuleDict['mlab'] = 'figureMlab'

    # list common arguments to pass
    figureConstructorSharedArgDict = {}
    figureConstructorSharedArgDict['msh'] = msh
    figureConstructorSharedArgDict['eqn'] = eqn

    # list figure specific arguments (can add figure-specific ones here)
    figureConstructorArgDict = {}
    mlabFigureArgDict = {}

    # add to general figure argument dictionary
    figureConstructorArgDict['mlab'] = mlabFigureArgDict

    figures = []
    for iB in (figureDict.keys()):

        figureType = figureDict.get(iB).get('type')

        thisFigureDict = figureDict.get(iB)
        # form arguments to pass to constructor for the i-th figure
        thisFigureArgs = {}
        thisFigureArgs.update(dict(figureDict.get(iB)))
       
        del thisFigureArgs['type']
        #print figureConstructorArgDict[figureTypeList[iB]]
        thisFigureArgs.update(figureConstructorArgDict.get(figureType))
        thisFigureArgs.update(figureConstructorSharedArgDict)
        thisFigureArgs['prefix'] = prefDir
        thisFigureArgs['name'] = iB

        if figureType in figureConstructorDict:
            f = __import__(figureModuleDict.get(figureType))
            figures.append(getattr(f,figureConstructorDict.get(figureType)) \
                (**thisFigureArgs))
            #figures.append(figureConstructorDict[figureType](**thisFigureArgs))
        else:
            sys.exit(' *** ERROR: Unsupported figure engine type! ***')

    #print figures
    return figures
##------------------------------------------------------------------------------
