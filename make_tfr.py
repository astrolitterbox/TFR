from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.fontsize'] = 10
#mpl.rcParams['font.size'] = 12
MLRatio = 1.7



def getHaloStellarMassRatio(stellarMassRange, alpha, beta, x0, y0, gamma):
    #ratio = log10y0*((np.divide(stellarMassRange, log10x0))**alpha)*(0.5+0.5*(np.divide(stellarMassRange, log10x0))**gamma)**((beta - alpha)/gamma)
    ratio1 = y0*(np.divide(stellarMassRange, x0))**alpha
    print ratio1
    print 'ratio1'
    ratio2 = 0.5+(0.5*(np.divide(stellarMassRange, x0)**gamma))**((beta-alpha)/gamma)
    print ratio2
    ratio = ratio1*ratio2
    return ratio


def getVvir(dynMass):
  #vVir = 1.5*(dynMass/(2.3*10**5))**(1./3) #1.5, from Mo, White, vdB, also McCarthy, Schaye
  vVir = (1/3)*np.log10((4.301*10**(-6))*dynMass)
  
  return vVir

def getAbsMag(MLRatio, Mstellar):
  print Mstellar, 'mstellar'
  L = (1/MLRatio)*(Mstellar)
  absMag = 4.76-2.5*np.log10(L)
  print absMag, 'absmag'
  return absMag


#coefficients from Dutton 2010b, Table 2, for late types, eq.3
alpha = -0.5
beta = 0.
x0 = 10**10.4
y0 = 10**1.61
gamma = 1.

alpha_lo = -0.45
y0_lo = 10**1.37

alpha_hi = -0.65
y0_hi = 10**1.89

#stellarMassRange = np.arange(0.7**(2)*10**9.3, 0.7**(2)*10**11.0, 10e8)
stellarMassRange = np.arange(10**9.3, 10**11.0, 10e8)

haloStellarMassRatio_lo = getHaloStellarMassRatio(stellarMassRange, alpha_lo, beta, x0, y0_lo, gamma)
haloStellarMassRatio_hi = getHaloStellarMassRatio(stellarMassRange, alpha_hi, beta, x0, y0_hi, gamma)
haloStellarMassRatio = getHaloStellarMassRatio(stellarMassRange, alpha, beta, x0, y0, gamma)

#repeating Figure 1:

fig = plt.figure(figsize=(8, 8))

ax = fig.add_subplot(221)
ax.plot(np.log10(stellarMassRange), np.log10(haloStellarMassRatio_lo), c='grey')
ax.plot(np.log10(stellarMassRange), np.log10(haloStellarMassRatio), c='r')
ax.plot(np.log10(stellarMassRange), np.log10(haloStellarMassRatio_hi), c='grey')
plt.fill_between(np.log10(stellarMassRange),np.log10(haloStellarMassRatio_lo), np.log10(haloStellarMassRatio_hi), color='grey', alpha='0.2')
plt.title('Stellar mass-halo ratio for late types')
plt.xlabel(r"$log(M_{*}), h^{-2}M_{\odot}$")
plt.ylabel(r"$log(M_{*})/log(M_{halo})$")

Mhalo_lo = stellarMassRange*haloStellarMassRatio_lo
Mhalo = stellarMassRange*haloStellarMassRatio
Mhalo_hi = stellarMassRange*haloStellarMassRatio_hi

stellarMassRange = stellarMassRange

ax = fig.add_subplot(222)
plt.title(' from Dutton 2010b')
ax.plot(np.log10(Mhalo_lo), stellarMassRange/Mhalo_lo, c='grey')
ax.plot(np.log10(Mhalo), stellarMassRange/Mhalo, c='r')
ax.plot(np.log10(Mhalo_hi), stellarMassRange/Mhalo_hi, c='grey')
ax.set_yscale('log')
plt.xlabel(r"$log(M_{halo}), h^{-1}M_{\odot}$")
plt.ylabel(r"$log(M_{*}),  h^{-2}M_{\odot}$")

dynMass_lo = Mhalo_lo + stellarMassRange
dynMass = Mhalo + stellarMassRange
dynMass_hi = Mhalo_hi + stellarMassRange

vVir_lo = getVvir(dynMass_lo)
vVir = getVvir(dynMass)
vVir_hi = getVvir(dynMass_hi)

vVir_ver = 2.185+0.281*(np.log10(stellarMassRange) - 10)
vVir_c = 2.145+0.259*(np.log10(stellarMassRange) - 10)

ax = fig.add_subplot(223)
plt.title('Stellar mass - circular velocity relation')
ax.plot(np.log10(stellarMassRange), vVir_lo, c='grey')
ax.plot(np.log10(stellarMassRange), vVir, c='r', label='Derived')
ax.plot(np.log10(stellarMassRange), vVir_hi, c='grey')

e = plt.plot(np.log10(stellarMassRange), vVir_c, c="g", label = "Courteau 1999")
e = plt.plot(np.log10(stellarMassRange), vVir_ver, c="b", label = "Verheijen 2001")
plt.fill_between(np.log10(stellarMassRange), vVir_lo, vVir_hi, color='grey', alpha='0.2')
plt.legend(loc="top left")
plt.ylabel(r"$log(v_{vir})$")
plt.xlabel(r"$log(M_{*}),  h^{-2}M_\odot$")


absMag = getAbsMag(MLRatio, stellarMassRange*0.7**-2)


ax = fig.add_subplot(224)
ax.plot(vVir, absMag)
plt.title('...Applying a M/L ratio:')
f = plt.plot(vVir, absMag, c="r", label = "Derived TFR")
g = plt.plot(vVir_c, absMag, c="g", label = "Courteau 1999")
f = plt.plot(vVir_lo, absMag, c="grey")
f = plt.plot(vVir_hi, absMag, c="grey")
g = plt.plot(vVir_ver, absMag, c="b", label = "Verheijen 2001")
plt.fill_betweenx(absMag, vVir_lo, vVir_hi, color='grey', alpha='0.2')
ax.axis([1.85, 2.6, -18.5, -23.0])
#e = plt.plot(vVir, -5.692*(vVir - 2.20) - 21.36, c="Green", label = "Pizagno 2007, slope:-5.7")

plt.legend(loc="best")
plt.xlabel(r"$log(v_{Vir})$")
plt.ylabel(r"$M_r, M_{*} (h = 0.7)$")

plt.tight_layout(pad=-0.4, w_pad=-2, h_pad=0.5)

plt.savefig('TFR', bbox_inches='tight')