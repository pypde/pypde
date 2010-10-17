#from copy import deepcopy

from probeForce import *

##------------------------------------------------------------------------------
def constructProbes(refState,msh,probeDict,prefDir):

    ###--------------------------------------------------------------------------
    ### set up lookup tables for class constructors
    ###--------------------------------------------------------------------------
    probeConstructorDict = {}
    probeConstructorDict['forces'] = ProbeForce

    # list common arguments to pass
    probeConstructorSharedArgDict = {}
    probeConstructorSharedArgDict['msh'] = msh
    

    # list probe specific arguments (can add probe-specific ones here)
    probeConstructorArgDict = {}
    forceProbeArgDict = {}
    forceProbeArgDict['refState'] = refState

    probeConstructorArgDict['forces'] = forceProbeArgDict

    # list probe specific arguments (can add equation specific ones here)
    #for iDict in probeConstructorArgDict.values():
    #    iDict.update(probeConstructorSharedArgDict)
 
    probes = []
    for iB in (probeDict.keys()):

        probeType = probeDict.get(iB).get('type')

        thisProbeDict = probeDict.get(iB)
        # form arguments to pass to constructor for the i-th probe
        thisProbeArgs = {}
        thisProbeArgs.update(dict(probeDict.get(iB)))
       
        del thisProbeArgs['type']
        #print probeConstructorArgDict[probeTypeList[iB]]
        thisProbeArgs.update(probeConstructorArgDict.get(probeType))
        thisProbeArgs.update(probeConstructorSharedArgDict)
        thisProbeArgs['prefix'] = prefDir
        thisProbeArgs['name'] = iB

        #print 'thisProbeArgs'
        #print thisProbeArgs
        probes.append(probeConstructorDict[probeType](**thisProbeArgs))

    #print probes
    return probes
##------------------------------------------------------------------------------
