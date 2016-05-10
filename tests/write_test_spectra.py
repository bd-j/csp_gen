import numpy as np
import fsps
import matplotlib.pyplot as pl
import time, sys, pickle
from itertools import product
from prospect.io.write_results import run_command

# need to also test tage=0 and SFH=0, 2, 3


cmd = 'cd $SPS_HOME; git rev-parse --abbrev-ref HEAD'
fsps_branch = (run_command(cmd)[1][0]).splitlines()[0]
sps = fsps.StellarPopulation(zcontinuous=1)


# want to test tau, tage, sf_trunc, sf_start, const, and fburst
default_sfh = {'tage': 10.0, 'tau': 3.0, 'const': 0.0, 'fburst': 0.0,
               'sf_start': 0.0, 'sf_trunc':0.0, 'sf_slope': 0.0, 'tburst':4.0}


def test_ssps(**pars):
    for k, v in pars.iteritems():
        sps.params[k] = v
    sps.params['sfh'] = 0
    blob = []
    ages = np.linspace(7, 10, 11)
    ages = np.insert(10**(ages-9), 0, 0)
    for tage in ages:
        w, s = sps.get_spectrum(tage=tage, peraa=True)
        blob.append(['tage', tage, s])
    return blob


def spardict(sps):
    return dict([(p, sps.params[p])
                 for p in sps.params.all_params])


def test_params_product(variations, sps=sps, **pars):
    parnames = np.sort(variations.keys())
    # Get all parameter combinations
    vpars = list(product(*[variations[n] for n in parnames]))
    ty = [type(np.array(variations[n])[0]) for n in parnames]
    parstruct = np.array(vpars, dtype=np.dtype(zip(parnames, ty)))

    # output blob
    all_blob = []
    spec = np.zeros([len(parstruct), len(sps.wavelengths)])
    mass = np.zeros(len(parstruct))
    
    # loop over parameter combinations
    for i,p in enumerate(parstruct):
        # Set any global parameters
        for k, v in pars.iteritems():
            sps.params[k] = v
        # set these specific parameters
        for n in parnames:
            sps.params[n] = p[n]
        # only proceed if we won't get an error for the params
        good = ((sps.params['tage'] - sps.params['sf_start'] > 0.1) &
                ((sps.params['fburst'] == 0) | (sps.params['const'] == 0))
                )
        if good:
            t = time.time()
            w, s = sps.get_spectrum(tage=sps.params['tage'], peraa=False)
            dt = time.time() - t
            m = sps.stellar_mass
        else:
            s = np.zeros(len(sps.wavelengths))
            m = 0
            dt = 0
        spec[i,:] = s
        mass[i] = m
        spars = spardict(sps)
        all_blob.append([spars])
    return [parstruct, spec, mass, all_blob]

            
def write_spectra(variations, sps=sps, sfh=1):
    pars = default_sfh.copy()
    pars['sfh'] = sfh
    filename = 'pickles/sfh{}_{}.pkl'.format(sfh, fsps_branch)
    blob = test_params_product(variations, **pars)
    blob.append([pars])
    with open(filename, 'wb') as f:
        pickle.dump(blob, f)

if __name__ == "__main__":

    # sfh = 0
    #for i, p in enumerate(parsets):
    #    filename = 'pickles/sfh0_pset{}_{}.pkl'.format(i, fsps_branch)
    #    blob = test_ssps(**p)
    #    with open(filename, 'wb') as f:
    #        pickle.dump(blob, f)

    #sys.exit()

    # sfh = 1
    variations = {'tage': np.linspace(0.2, 10, 4),
                  'tau': 10**np.linspace(-1, 2, 4),
                  'const': np.linspace(0.0, 1.0, 3),
                  'fburst': [0, 0.5],
                  'sf_start': default_sfh['tage'] * np.linspace(0.0, 0.9, 4),
                  'dust2': [0.0, 1.0],
                  'add_neb_emission': [0, 1],
                  'sf_trunc': np.linspace(0, 10, 3)
                  }

    t = time.time()
    write_spectra(variations, sfh=1)
    dt = time.time() - t
    
    # sfh = 4
    variations = {'tage': np.linspace(0.2, 10, 4),
                  'tau': 10**np.linspace(-1, 2, 4),
                  'const': np.linspace(0.0, 1.0, 3),
                  'fburst': [0, 0.5],
                  'sf_start': default_sfh['tage'] * np.linspace(0.0, 0.9, 4),
                  'dust2': [0.0, 1.0],
                  'add_neb_emission': [0, 1],
                  'sf_trunc': np.linspace(0, 10, 3)
                  }

    write_spectra(variations, sfh=4)

    # sfh = 5
    variations = {'tage': np.linspace(0.2, 10, 10),
                  'tau': np.linspace(0.1, 100, 3),
                  'sf_start': default_sfh['tage'] * np.linspace(0.0, 0.9, 4),
                  'sf_trunc': np.linspace(0, 10, 5),
                  'sf_slope': [-0.5, 0., 0.5]
                  }

    write_spectra(variations, sfh=5)
