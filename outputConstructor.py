##------------------------------------------------------------------------------
def constructOutput(eqn,msh,outputDict,prefDir):

    ###--------------------------------------------------------------------------
    ### set up lookup tables for class constructors
    ###--------------------------------------------------------------------------
    outputConstructorDict = {}
    outputConstructorDict['plain'] = 'OutputPlain'

    outputModuleDict = {}
    outputModuleDict['plain'] = 'outputPlain'

    # list common arguments to pass
    outputConstructorSharedArgDict = {}
    outputConstructorSharedArgDict['msh'] = msh
    outputConstructorSharedArgDict['eqn'] = eqn

    # list output specific arguments (can add output-specific ones here)
    outputConstructorArgDict = {}
    plainOutputArgDict = {}

    # add to general output argument dictionary
    outputConstructorArgDict['plain'] = plainOutputArgDict

    outputArgs = {}
    outputArgs.update(outputDict)
  
    outputType = outputDict.get('format')
    del outputArgs['format']
    del outputArgs['writeInterval']
    del outputArgs['figureInterval']
    del outputArgs['makeFigures']
    del outputArgs['figures']

    outputArgs.update(outputConstructorArgDict.get(outputType))
    outputArgs.update(outputConstructorSharedArgDict)
    outputArgs['prefix'] = prefDir

    if outputType in outputConstructorDict:
        f = __import__(outputModuleDict.get(outputType))
        output = getattr(f,outputConstructorDict.get(outputType)) \
            (**outputArgs)
    else:
        sys.exit(' *** ERROR: Unsupported output engine type! ***')

    #print outputs
    return output
##------------------------------------------------------------------------------
