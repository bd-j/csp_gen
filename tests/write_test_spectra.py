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
               'sf_start': 0.0, 'sf_trunc':0.0, 'sf_slope': 0.0, 'tburst':4.0,
               'add_neb_emission': False, 'dust2': 0.0}


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

    # output blobs
    all_blob, spec, mass, sfr, goodpars = [], [], [], [], []
    
    # loop over parameter combinations
    for i,thisp in enumerate(vpars):
        p = dict(zip(parnames, thisp))
        # Set any global parameters
        for k, v in pars.iteritems():
            sps.params[k] = v
        # set these specific parameters
        for n in parnames:
            sps.params[n] = p[n]
        # only proceed if we won't get an error for the params
        good = (
                # dies otherwise
                ((sps.params['tage'] - sps.params['sf_start'] > 0.1) | (sps.params['tage'] == 0)) &
                # both nonzero behaves differently by design
                ((sps.params['fburst'] == 0) | (sps.params['const'] == 0)) &
                # sf_start > sf_trunc behaves differently by design
                ((sps.params['sf_start'] <= sps.params['sf_trunc']))
                )
        if good:
            
            t = time.time()
            w, s = sps.get_spectrum(tage=sps.params['tage'], peraa=False)
            dt = time.time() - t
            mass.append(sps.stellar_mass)
            sfr.append(sps.sfr)
            spec.append(s)
            spars = spardict(sps)
            all_blob.append([spars])
            goodpars.append(thisp)

    print(len(goodpars))
    parstruct = np.array(goodpars, dtype=np.dtype(zip(parnames, ty)))
    return [parstruct, np.array(spec), np.array(mass), np.array(sfr), all_blob]

            
def write_spectra(variations, sps=sps, sfh=1, root='pickles/'):
    ages = sps.ssp_ages
    pars = default_sfh.copy()
    pars['sfh'] = sfh
    filename = '{}sfh{}_{}.pkl'.format(root, sfh, fsps_branch)
    blob = test_params_product(variations, **pars)
    blob.extend([pars, ages])
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


    #tage = 0
    variations = {'tage': [0.0],
                  'tau': 10**np.linspace(-1, 2, 4),
                  'const': [0.0, 0.5],
                  'fburst': [0, 0.5],
                  'sf_start': [0.0, 5],
                  'sf_trunc': np.linspace(0, 10, 4)
                  }
    write_spectra(variations, sfh=4, root='pickles/tage0_')

    variations = {'tage': [0.0],
                  'tau': 10**np.linspace(-1, 2, 4),
                  'sf_start': [0, 5],
                  'sf_trunc': np.linspace(0, 10, 4),
                  'sf_slope': [-0.5, 0.0, 0.5]
                  }
    t = time.time()
    write_spectra(variations, sfh=5, root='pickles/tage0_')
    dt = time.time() - t
    print(dt)
    
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

    write_spectra(variations, sfh=1)
    
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
                  'tau': 10**np.linspace(-1, 2, 4),
                  'sf_start': default_sfh['tage'] * np.linspace(0.0, 0.9, 4),
                  'sf_trunc': np.linspace(0, 10, 6),
                  'sf_slope': [-0.5, 0., 0.5],
                  'add_neb_emission': [0, 1],
                  }

    t = time.time()
    write_spectra(variations, sfh=5)
    dt = time.time() - t
