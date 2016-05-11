import numpy as np
import fsps
import matplotlib.pyplot as pl
import time

sps = fsps.StellarPopulation()

sps.params['sfh'] = 1
sps.params['tau'] = 1.0
sps.params['sf_trunc'] = 2.0
sps.params['dust2'] = 0.0
sps.params['add_neb_emission'] = True
sps.params['add_neb_continuum'] = False
sps.params['add_dust_emission'] = False
w, s = sps.get_spectrum(tage=0, peraa=True)
pl.plot(w, s[-1,:])
pl.xscale('log')
pl.yscale('log')
pl.show()

sps.params['const'] = 0.5
for tage in np.linspace(10, 0.11, 10):
    t =time.time()
    w, s = sps.get_spectrum(tage=tage, peraa=True)
    dt = time.time() - t
    print(tage, dt)
print('doing full')
for tau in np.linspace(1, 2, 10):
    sps.params['tau'] = tau
    t =time.time()
    w, s = sps.get_spectrum(tage=0, peraa=True)
    dt = time.time() - t
    print(dt)

for d in np.linspace(0, 1, 10):
    sps.params['sfh'] = 0
    sps.params['dust2'] = d
    t =time.time()
    w, s = sps.get_spectrum(tage=1, peraa=True)
    dt = time.time() - t
    print(dt)
