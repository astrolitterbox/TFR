#A pet project of mine -- taking the parameterised form of M*-M_halo relation and converting it to luminosity-rotation velocity relation
#i.e. looking for the origin of the Tully-Fisher relation
#Relevant citations: Vogelsberger et al. http://arxiv.org/pdf/1305.2913v1.pdf, p.21
#Behroozi 2013, Conroy 2006, Mo, White, van den Bosch 2000

#McGaugh 2000. de Jong (1996) -- M*/L ratio: 

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def getStellarMass(Mhalo, Mratio):
   Mstellar = np.log10(np.multiply(Mhalo, Mratio))
   return Mstellar

def getAbsMag(MLRatio, Mstellar):
  L = MLRatio*10**(Mstellar)
  absMag = -2.5*np.log10(L)
  return absMag

def getVvir(galaxies, dynMass):
  vVir = 1.5*(dynMass/(2.3*10**5))**(1./3) #1.5, from Mo, White, vdB
  vVir = vVir[galaxies]
  return vVir


MLRatio = 1/1.7 # I band

data = np.genfromtxt("behroozi.csv",  usecols=(0, 1, 2, 3))
Mhalo = 10**data[:, 0]
Mratio = 10**data[:, 1]
Mratio_lo = 10**(data[:, 1] - data[:, 3])

Mratio_hi = 10**(data[:, 1] + data[:, 2])
#print Mratio
Mstellar_hi = getStellarMass(Mhalo, Mratio_hi)
Mstellar_lo = getStellarMass(Mhalo, Mratio_lo)
Mstellar = getStellarMass(Mhalo, Mratio)
#Mhalo = np.log10(Mhalo)


fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
plt.plot(Mhalo, Mratio)
plt.plot(Mhalo, Mratio_lo)
plt.plot(Mhalo, Mratio_hi)
plt.savefig("M_halo_M_s")


absMag = getAbsMag(MLRatio, Mstellar)
galaxies = np.where(absMag > -24)
absMag = absMag[galaxies]
absMag_lo = getAbsMag(MLRatio, Mstellar_lo)[galaxies]
absMag_hi = getAbsMag(MLRatio, Mstellar_hi)[galaxies]

dynMass = Mhalo + 10**Mstellar
dynMass_lo = Mhalo + 10**Mstellar_lo
dynMass_hi = Mhalo + 10**Mstellar_hi

print dynMass

#M = 2.33*10^5*vVir^3
#vVir = np.log10(dynMass**(0.3))
vVir = getVvir(galaxies, dynMass)
vVir_lo = getVvir(galaxies, dynMass_lo)
vVir_hi = getVvir(galaxies, dynMass_hi)

fig = plt.figure()
ax = fig.add_subplot(111)
e = plt.plot(np.log10(vVir), -5.55*(np.log10(vVir) - 2.20) - 21.36, c="k", label = "Pizagno 2007")

ax.plot(np.log10(vVir), absMag)
ax.plot(np.log10(vVir_lo), absMag_lo)
ax.plot(np.log10(vVir_hi), absMag_hi)
plt.ylim(plt.ylim()[::-1])
plt.savefig("TF")

