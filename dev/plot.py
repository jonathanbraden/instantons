#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from palettable import *
#import matplotlib.cm as cm
#import matplotlib.colors as colors

files = [ 'instanton_{0:}.dat'.format(l) for l in ['0.50001', '0.5001', '0.501', '0.51', '0.52', '0.6', '0.8', '1', '2', '4', '6', '8', '50' ] ]
files_2 = [ 'i_{0:}.dat'.format(l) for l in [ '0.55', '0.58', '0.6', '0.7', '0.8', '1', '5', '10', '50' ] ]

actionOrders = [ 10, 25, 50, 100, 200 ]
actionFiles = [ 'actions-o{0:}.dat'.format(str(o)) for o in actionOrders ]

def plotProfiles(fNames):
    cm = plt.get_cmap('OrRd')
    colors = cm( np.linspace(0.2,1,len(files)) )
    rInd = 0; fInd = 2 # Index for radius and field values
    
    for c,f in zip(colors,fNames):
        a = np.genfromtxt(f)
        plt.plot(a[:,rInd],a[:,fInd],color=c)
    plt.xlabel(r'$mr_{\rm E}$')
    plt.ylabel(r'$\phi$')
    plt.xlim(0,25)

def plotAction(file):
    a = np.genfromtxt(file)

    ind = [1,3,4,5,7,6]
    labels = [ r'$S_{\rm E}$', r'$\frac{2}{D}S_{\rm KE}$', r'$-\frac{D-2}{D}S_{\rm PE}$', r'$\frac{2}{D}R_{\rm tw}^d S_{1}$', r"$\frac{1}{D}S_{\phi V'}$", r'$S_{V_{\rm tw}}$' ]
    dim=4
    norms = [ 1., 2./dim, -(dim-2.)/dim, 0.5, 0.25, 1. ]
    for j,(i,l) in enumerate(zip(ind,labels)):
        plt.plot(a[:,0]-0.5,norms[j]*np.abs(a[:,i]),label=l)
    plt.xlabel(r'$\delta-\frac{1}{2}$'); plt.ylabel(r'$S / \Omega_D$')
    plt.xscale('log'); plt.yscale('log')
    #plt.legend(loc='upper left',bbox_to_anchor=(0,0,1,1),bbox_transform=plt.gca().transAxes)
    return

def compareActions(files,orders,ind,norm=1.):
    for o,f in zip(orders,files):
        a = np.genfromtxt(f)
        plt.plot(a[:,0]-0.5,norm*a[:,ind],label=r'${0:}$'.format(o))
    plt.xlabel(r'$\delta-\frac{1}{2}$')
    plt.legend(loc='upper left', bbox_to_anchor=(0,0,1,1),bbox_transform=plt.gca().transAxes,title=r'$N_{\rm max}$')
    plt.xscale('log'); plt.yscale('log')
    return
    
def plotSurfaceTension(file):
    a = np.genfromtxt(file)

    ind = 2
    plt.plot(a[:,0]-0.5,a[:,ind])
    plt.xlabel(r'$\delta-\frac{1}{2}$'); plt.ylabel(r'$\sigma$')
    plt.xscale('log'); plt.yscale('log')
    return

if __name__=="__main__":
    plotProfiles(files)
    plt.xlim(0.,40.)
    plt.savefig('instanton-profiles.pdf'); plt.savefig('instanton-profiles.png')
    
