#import pylab
#import operator
import sys
import os
from enthought.mayavi import mlab
from numpy import *
#from numpy.linalg import *
#from pysparse import spmatrix
#from pysparse import itsolvers
from pyparsing import *

from caseManager import *
from eulerEquations import *
from schmSpaceConstructor import *
from solverConstructor import *
from outputConstructor import *
from figureConstructor import *
from probeConstructor import *

##==============================================================================
## Function defn
##==============================================================================
#-------------------------------------------------------------------------------
#def getDeltaT(cfl,tModeLocal,w,eqn,msh):
# 
#    dt = cfl*msh.mass / eqn.spectralRadius(w,msh)
#
#    if not(tModeLocal):
#        dt = min(dt)*ones(dt.shape)
#        dt = min(dt)*ones(msh.nP)
#    return dt
#-------------------------------------------------------------------------------
#def outputFig(msh,f,s,it,sName,sDir):
#    mlab.clf()
#    mlab.points3d(msh.x[:,0], msh.x[:,1], f, s, 
#                  scale_factor=.05, scale_mode='none')
#    mlab.outline(extent=[msh.BB[0,0],msh.BB[1,0],msh.BB[0,1],msh.BB[1,1],0,0])
#    #mlab.outline(extent=[-.2, 1.2, -.2, 1.2,0,0])
#    mlab.view(0,0,5,[0.5,0,0])
#    #mlab.view(0,0,50,[0.5,0,0])
#    #mlab.view(0,0,150,[0.5,0,0])
#    mlab.colorbar()
#    fSol =  ''.join([sDir,'/',sName,'%04d.png'% it])
#    mlab.savefig(fSol)
##-------------------------------------------------------------------------------
#def outputSol(msh,w,it,sDir):
#    f = open(''.join([sDir,'/sol%04d.txt' % it]),'w')
#    for i in range(msh.nP):
#        (w[:,i]).tofile(f,format='%20.12e', sep=' ')
#        f.write('\n')
#    f.close
##-------------------------------------------------------------------------------

##==============================================================================
## Main program
##==============================================================================

print "Starting simulation"
# get command line arguement for directory and so on

argList = sys.argv
nArgs = len(argList)
iargc = 1 # start from the second argument to skip the executable name
prefDir = os.getcwd()
fileIn = 'problemDefinition'
# NEED TO CHECK DIRECTORY CONTAIN VALID FILES.
hasPrefix = False

while iargc < len(argList):
    arg = argList[iargc]
    if arg == '-dir':
       if not hasPrefix:
           iargc += 1
           prefDir = argList[iargc]
           hasPrefix == True
    elif arg == '-file':
           iargc += 1
           fileIn = argList[iargc]
    else:
       sys.exit('Invalid command-line option.' \
                'The only one supported is prefix directory')
    iargc +=1

print ' ... Prefix is:', prefDir

# set variables
solDir = ''.join([prefDir,'/solutions'])
if not(os.path.isdir(solDir)):
    os.mkdir(solDir)
#if 'cleanDirectory' in case.outputDict:	# DO NOTHING EVEN IF ASKED TO WIPE FILES -- TOO DANGEROUS
    #if case.outputDict['cleanDirectory']:
    #    os.system(''.join(['rm ',solDir,'/*']))

##------------------------------------------------------------------------------
## Read input file
##------------------------------------------------------------------------------
# in the future, do this for each subdomain

case = CaseManager(prefDir,fileIn)
msh = case.msh

# iteration control
nt   = case.timeDict['timeSteps']
del case.timeDict['timeSteps']

nWrite = case.outputDict['writeInterval']
nFig = case.outputDict['figureInterval']
makeFigs = case.outputDict['makeFigures']

##------------------------------------------------------------------------------
## initialize variables
##------------------------------------------------------------------------------
print "------------------------------------------------------------------------"
print " Case parameters"
print "------------------------------------------------------------------------"

# build equation model
interfaceSolverDict = case.spaceDict.get('interfaceSolvers')

# this is an Euler solver...so just initiate Euler equation directly
eqn = EulerEquations(msh,
                     case.modelDict,
                     case.condsDict.get('reference'),
                     case.bcDict,
                     interfaceSolverDict)
# if equation is chosen based on input file:
#eqn = constructEqn(msh,**case.modelDict)

# build spatial scheme to solve the equations
spaceSchm = constructSchmSpace(eqn,msh,case.spaceDict)

# build temporal scheme to solve the equations
solver = constructSolver(case.timeDict)

# output
output = constructOutput(eqn,msh,case.outputDict,prefDir)

# figures
figuresList = []
if case.outputDict['makeFigures']:
	figuresList = constructFigures(eqn,msh,case.outputDict.get('figures'),prefDir)

# force probes
probesList = constructProbes(eqn.refState(),msh,case.probesDict,prefDir)

##------------------------------------------------------------------------------


##------------------------------------------------------------------------------
# print parameters
##------------------------------------------------------------------------------
print " ... Time iteration control"
print "        number of time steps:   ", nt
print "        write interval:         ", nWrite
print

##------------------------------------------------------------------------------
## initialize solution in domain
##------------------------------------------------------------------------------
# uniform flow from reference state
eqn.initSolution()

##------------------------------------------------------------------------------
## alias solutions
##------------------------------------------------------------------------------
# write initial solution
output.save(0)
if makeFigs:
    for fig in figuresList:
	   fig.save(0)
##------------------------------------------------------------------------------

fnameRes = ''.join([prefDir,'/','residual.txt'])
fileRes = open(fnameRes,'w')
fileRes.write('# it maxRes(mass) maxRes(x-mom) maxRes(y-mom) maxRes(energy)\n')
fileRes.close()

it = 0
converged = False

while((it < nt) and (not(converged))):
#for it in range(nt):

    converged, maxRes = solver.solve(eqn,spaceSchm)

    print it, maxRes
    sys.stdout.flush()

    fileRes = open(fnameRes,'a')
    fileRes.write('%d' % (it) )
    maxRes.tofile(fileRes,format='%24.12e', sep=' ')
    fileRes.write('\n')

    # print solution at selected intervals
    if ((it+1) % nWrite) == 0:
        ## calculate Mach
        #localM = eqn.mach()
        #outputSol(msh,vstack((eqn.solution(),
        #                      eqn.extraVariables(),
        #                      localM)),
        #          it+1,solDir)
        output.save(it+1)
        for probes in probesList:   
        	probes.writePressure(it+1,eqn.allVariables(),msh)
    if ((it+1) % nFig) == 0:
        if makeFigs:
            for fig in figuresList:
        	   fig.save(it+1)
        #if makeFigs:
        #    localM = eqn.mach()
        #    outputFig(msh,localM,localM,it+1,'M',solDir)
    # probe stuff at every iteration
    for probes in probesList:   
	probes.writeProbeValues(it+1,eqn.allVariables(),msh)
    
    it = it+1 
    #...........................................................................
    # end of time iterations 
    #...........................................................................

# final output    
# calculate Mach
if ((it) % nWrite) != 0:
    #localM = eqn.mach()
    #outputSol(msh,vstack((eqn.allVariables(),localM)),it,solDir)
    output.save(it)
    for probes in probesList:   
        probes.writePressure(it,eqn.allVariables(),msh)

if ((it) % nFig) != 0:
    for fig in figuresList:
        fig.save(it)
    #if makeFigs:
    #   outputFig(msh,localM,localM,it,'M',solDir)
## end of program
