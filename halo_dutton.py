
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

likkle_h = 0.7

#From Dutton 2010b:

def getHaloStellarMassRatio(stellarMassRange, alpha, beta, x0, y0, gamma):
    #ratio = log10y0*((np.divide(stellarMassRange, log10x0))**alpha)*(0.5+0.5*(np.divide(stellarMassRange, log10x0))**gamma)**((beta - alpha)/gamma)
    ratio1 = y0*(np.divide(stellarMassRange, x0))**alpha
    print ratio1
    print 'ratio1'
    ratio2 = 0.5+(0.5*(np.divide(stellarMassRange, x0)**gamma))**((beta-alpha)/gamma)
    print ratio2
    ratio = ratio1*ratio2
    return ratio
    
def getStellarMass(Mhalo, Mratio):
   Mstellar = np.log10(np.multiply(Mhalo, Mratio))
   return Mstellar


def getAbsMag(MLRatio, Mstellar):
  print Mstellar, 'mstellar'
  L = (1/MLRatio)*(10**(Mstellar))
  absMag = 4.77-2.5*np.log10(L)
  print absMag, 'absmag'
  return absMag

def getVvir(dynMass):
  #vVir = 1*(dynMass/(2.3*10**5))**(1./3) #1.5, from Mo, White, vdB, also McCarthy, Schaye
  vVir = 10**(1/3)*np.log10((4.301*10**(-6))*dynMass)
  return vVir


stellarMassRange = np.arange(0.7**(2)*10**9.3, 0.7**(2)*10**11.0, 10e8)

#Mass to light ratio:
MLRatio = 1.7

#coefficients from Dutton 2010b, Table 2, for late types
alpha = -0.5
beta = 0.
x0 = 10**10.4
y0 = 10**1.61
gamma = 1.

alpha_lo = -0.45
y0_lo = 10**1.37

alpha_hi = -0.65
y0_hi = 10**1.89


haloStellarMassRatio_lo = getHaloStellarMassRatio(stellarMassRange, alpha_lo, beta, x0, y0_lo, gamma)
haloStellarMassRatio_hi = getHaloStellarMassRatio(stellarMassRange, alpha_hi, beta, x0, y0_hi, gamma)
haloStellarMassRatio = getHaloStellarMassRatio(stellarMassRange, alpha, beta, x0, y0, gamma)

#repeating Figure 1:

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(221)
ax.plot(np.log10(stellarMassRange), np.log10(haloStellarMassRatio_lo))
ax.plot(np.log10(stellarMassRange), np.log10(haloStellarMassRatio_hi))
ax.plot(np.log10(stellarMassRange), np.log10(haloStellarMassRatio))
plt.title('Stellar mass-halo ratio, Dutton')
plt.xlabel("log(M_star)")
plt.ylabel("log(M_star)/log(M_halo)")

print np.log10(haloStellarMassRatio), 'ratio'
print np.log10(stellarMassRange), 'stellar'

#M_stellar - M_halo relation:
haloMass = np.multiply(haloStellarMassRatio,stellarMassRange)
haloMass_lo = np.multiply(haloStellarMassRatio_lo,stellarMassRange)
haloMass_hi = np.multiply(haloStellarMassRatio_hi,stellarMassRange)

print haloMass


#From Behroozi 2012:
data = np.genfromtxt("behroozi.csv",  usecols=(0, 1, 2, 3))
Mhalo = 10**data[:, 0]
Mratio = 10**data[:, 1]
Mstellar = getStellarMass(Mhalo, Mratio)

ax = fig.add_subplot(222)
ax.plot(np.log10(haloMass), np.log10(stellarMassRange), label='Dutton', c='r')
ax.plot(np.log10(haloMass_lo), np.log10(stellarMassRange), c='g')
ax.plot(np.log10(haloMass_hi), np.log10(stellarMassRange), c='g')
ax.plot(np.log10(Mhalo), Mstellar, c='b', label='Behroozi')
plt.xlabel("log(M_halo)")
plt.ylabel("log(M_star)")

plt.legend()
#******************************************************8
#Trying to obtain the TFR


#cropping Behroozi data:
absMag_B = getAbsMag(MLRatio, Mstellar)
galaxies = np.where(absMag_B <> -24)
absMag_B = absMag_B[galaxies]
Mhalo = Mhalo[galaxies]
Mstellar = Mstellar[galaxies]
#Luminosities for Dutton
absMag_D = getAbsMag(MLRatio, np.log10(stellarMassRange))


#virial velocities for Behroozi 2012 data
dynMass = Mhalo + 10**Mstellar
vVir_B = getVvir(dynMass)

#virial velocities for Dutton data
vVir_D = getVvir(haloMass)
vVir_D_lo = getVvir(haloMass_lo)
vVir_D_hi = getVvir(haloMass_hi)


ax = fig.add_subplot(223)
ax.plot(np.log10(vVir_B), absMag_B, label='Behroozi', c='b')
ax.plot(np.log10(vVir_D), absMag_D, c='red', label='Dutton')
ax.plot(np.log10(vVir_D_hi), absMag_D, c='g')
ax.plot(np.log10(vVir_D_lo), absMag_D, c='g')

e = plt.plot(np.log10(vVir_B), -5.55*(np.log10(vVir_B) - 2.20) - 21.36, c="g", label = "Pizagno 2007")
e = plt.plot(np.log10(vVir_B), 2.143 + 0.281*(np.log10((10**Mstellar)*(1/MLRatio)/10**10)), label = "Pizagno 2007")
plt.legend()
plt.xlabel("log(1.2*v_vir)")
plt.ylabel("M_r")


plt.savefig('TF')
