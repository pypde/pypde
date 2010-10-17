from numpy import *

class RefStateFreeStream:
    #def __init__(self,gasModel,**kwarg):
        #self.alpha = kwarg['alpha']
        #self.mach = kwarg['mach']

    def __init__(self,gasModel, alpha=0., mach=0.775, cl=None):

        #print alpha
        #print mach
        self.__gma = gasModel.gma()
        self.__alpha = alpha
        self.__ralpha = alpha*pi/180.
        self.__mach = mach
        self.__cl = cl

        self.__p = 1.
        self.__rho = 1.
        #self.u0 = mach*c0
        #self.v0 = 0.
        #self.ei0 = p0/((gma-1)*rho0)
        #self.h0 = gma*ei0 + .5*(u0*u0 + v0*v0)

    def gma(self):
        return self.__gma

    def mach(self):
        return self.__mach

    def alpha(self):
        return self.__alpha

    def c(self):
        return sqrt(self.__gma*self.__p/self.__rho)

    def p(self):
        return self.__p

    def rho(self):
        return self.__rho

    def u(self):
        return self.__mach * self.c() * cos(self.__ralpha)

    def v(self):
        return self.__mach * self.c() * sin(self.__ralpha)

    def ei(self):
        return self.__p/((self.__gma-1.)*self.__rho)

    def h(self):
        return self.__gma* self.ei() + .5*((self.u())**2 + (self.v())**2)

    def s(self):
        return ((self.__rho**(self.__gma))/ self.__p)

    def velMag(self):
        return sqrt((self.u())**2. + (self.v())**2.)

    def setAlpha(self,alpha):
        self.__ralpha = alpha*(pi/180.)
        self.__alpha = alpha

    def setAlphaRadians(self,ralpha):
        self.__alpha = ralpha*(180./pi)
        self.__ralpha = ralpha
