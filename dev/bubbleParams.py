#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, brentq
from scipy.integrate import quad, odeint

class BubbleParams(object):
    def __init__(self,pot,dpot,phis):
        self.phiFalse = fsolve(dpot,phis[0])[0]
        self.phiTrue  = fsolve(dpot,phis[1])[0]
        Vout = lambda x : pot(x) - pot(self.phiFalse)
        self.phiOut = brentq(Vout,np.min([self.phiFalse,self.phiTrue])+0.1,np.max([self.phiFalse,self.phiTrue])-0.1)
        
        self.sigma = quad(pot,self.phiFalse,self.phiOut)[0]
        self.dRho = pot(self.phiFalse) - pot(self.phiTrue)
        self.rInit = self.sigma / self.dRho # R_init / dimension

        self.V = pot
        self.Vp = dpot

    def __str__(self):
        out = "False vacuum is {0:} and true vacuum is {1:0}".format(self.phiFalse,self.phiTrue)
        out +="\nTunnel out location is {0:}".format(self.phiOut)
        out += "\nTension is {0:} and delta rho is {1:}".format(self.sigma,self.dRho)
        out += "\nInitial bubble radius is (d-1){0:}".format(self.rInit)
        return out

    def plot(self):
        dPhi = np.abs(self.phiTrue-self.phiFalse)
        # This assumes phiTrue > phiFalse
        xvals = np.linspace(self.phiFalse-0.5*dPhi,self.phiTrue+0.5*dPhi,101)
        plt.plot(xvals,self.V(xvals))
        plt.plot(self.phiFalse,self.V(self.phiFalse),'go')
        plt.plot(self.phiTrue,self.V(self.phiTrue),'bo')
        plt.plot(self.phiOut,self.V(self.phiOut),'ro')
        plt.xlabel(r'$\phi$'); plt.ylabel(r'$V(\phi)$')
        plt.show()
        
    def testCalc(self):
        print "Derivative a false vacuum is ",self.Vp(self.phiFalse)
        print "Derivative a true vacuum is ",self.Vp(self.phiTrue)
        print "Potential diff at phi out is ",self.V(self.phiOut)-self.V(self.phiFalse)

        # This is unnecessary with the bracketed root finder above
        if (self.phiTrue > self.phiFalse):
            if (self.phiOut > self.phiTrue):
                print "The tunnel out point is past the true vacuum"
        else:
            if (self.phiOut < self.phiTrue):
                print "The tunnel out point is past the true vacuum"

        if (abs(self.phiTrue-self.phiFalse) < 0.01):
            print "Warning, false and true vacua are less that 0.01 apart"

    def profile(self):
        """
        Compute an approximate thin-wall profile
        """
        sol = odeint()
        
if __name__=="__main__":
    V = lambda x,e : 0.25*(x**2-1)**2 - e*x
    Vp = lambda x,e : x*(x**2-1) - e
    phi0 = [-1,1]
    
    eps = 0.01
    Bubble = BubbleParams(lambda x : V(x,eps), lambda x : Vp(x,eps), phi0)

    V = lambda x,e : -np.cos(x) + e*np.sin(x)**2
    Vp = lambda x,e : np.sin(x) + e*np.sin(2.*x)
    phi0 = [np.pi,2*np.pi]
    eps = 1.
    Bubble = BubbleParams(lambda x : V(x,eps), lambda x : Vp(x,eps), phi0)
