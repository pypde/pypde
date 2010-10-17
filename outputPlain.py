import os
from numpy import *

class OutputPlain(object):
    ##--------------------------------------------------------------------------
    def __init__(self, prefix,
                       msh=None,
                       eqn=None,
                       variables=None):

        self.__msh = msh
        self.__eqn = eqn
        self.__varList = variables
        self.__nVar = len(variables)

        # create directory
        self.__writeDir = ''.join([prefix,'/solutions/plain'])
        #try:
        #    os.makedirs(self.__writeDir)
        #except OSError as exc: # Python >2.5
        #    if exc.errno == errno.EEXIST:
        #        pass
        #    else: raise

        if not(os.path.isdir(self.__writeDir)):
            os.makedirs(self.__writeDir)
       
    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    def save(self,it):
    #def outputSol(self,msh,w,it,sDir):
        f = open(''.join([self.__writeDir,'/sol%04d.txt' % it]),'w')
        w = empty([self.__nVar,self.__msh.nP])
        for i in range(self.__nVar):
            w[i,:] = (self.__eqn).get(self.__varList[i])()
        
        for i in range(self.__msh.nP):
            w[:,i].tofile(f,format='%20.12e', sep=' ')
            f.write('\n')
        f.close

