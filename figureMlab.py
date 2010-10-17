import os
from numpy import *
from enthought.mayavi import mlab

class FigureMlab(object):
    ##--------------------------------------------------------------------------
    def __init__(self, prefix,
                       msh=None, 
                       eqn=None, 
                       name='field',
                       variable=None,
                       glyphScaleFactor=0.1, glyphScaleMode='none',
                       viewZoom=10., viewCenter=[0.,0.,0.]):

        self.__msh = msh
        self.__eqn = eqn
        self.__var = variable
        self.__name = name
        self.__mode = glyphScaleMode
        self.__scale = glyphScaleFactor
        self.__zoom = viewZoom
        self.__xView = list(viewCenter)

        # create directory
        self.__figDir = ''.join([prefix,'/solutions/figs/',name])
        #try:
        #    os.makedirs(self.__figDir)
        #except OSError as exc: # Python >2.5
        #    if exc.errno == errno.EEXIST:
        #        pass
        #    else: raise

        if not(os.path.isdir(self.__figDir)):
            os.makedirs(self.__figDir)
       
    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    def save(self,it):

        msh = self.__msh
        w = (self.__eqn).get(self.__var)()

        mlab.clf()
        mlab.points3d(msh.x[:,0], msh.x[:,1], zeros(w.shape), w,
                      scale_factor=self.__scale, scale_mode=self.__mode)
        mlab.outline(extent=[msh.BB[0,0],msh.BB[1,0],msh.BB[0,1],msh.BB[1,1],
                     0,0])
        mlab.view(0,0,self.__zoom,self.__xView)
        mlab.colorbar()
        fSol =  ''.join([self.__figDir,'/',self.__name,'%04d.png'% it])
        #print fSol
        mlab.savefig(fSol)
    ##--------------------------------------------------------------------------
