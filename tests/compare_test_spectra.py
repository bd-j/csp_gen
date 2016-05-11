import pickle, sys
import numpy as np
import matplotlib.pyplot as pl
import fsps
sps = fsps.StellarPopulation(zcontinuous=1)

sfh = 5
root = 'pickles/'

def main(sfh=1, root='pickles/'):
    pars, specs, masses, sfrs = read_pickles(sfh=sfh, root=root)    
    #missing = [np.all(branch_spec == 0., axis=-1) for branch_spec in spec]
    #either_missing = missing[0] | missing[1]
    #both_missing = missing[0] & missing[1]
    return plot_ratios(pars, specs, masses, sfh=sfh)
    
        
def read_pickles(sfh, branches=['localvars', 'new_compsp'], root='pickles/'):
    spec, parstruct, masses, sfrs = [], [], [], []
    print(branches)
    for branch in branches:
        print(branch)
        filename = '{}sfh{}_{}.pkl'.format(root, sfh, branch)
        with open(filename, 'rb') as f:
            blob = pickle.load(f)
        print(len(blob))
        try:
            branch_pars, branch_spec, m, sf, pardicts, default_sfh = blob
        except(ValueError):
            branch_pars, branch_spec, m, pardicts, default_sfh = blob
            sf = None
        spec.append(branch_spec)
        parstruct.append(branch_pars)
        masses.append(m)
        sfrs.append(sf)
    return parstruct, spec, masses, sfrs


def plot_sfh():
    pass


def plot_ratios(pars, specs, masses, sfh=4, sps=sps,
                ishow=[100, 500, 700], wlo=1e3, whi=2e4):

    wave = sps.wavelengths
    imin, imax = np.argmin(np.abs(wave - wlo)), np.argmin(np.abs(wave - whi)) 
    sratio = specs[0]/ specs[1]
    mratio = masses[0] / masses[1]
    maxratio = np.max(np.abs(np.log(sratio[:, imin:imax])), axis=-1)
    maxratioloc = np.argmax(np.abs(np.log(sratio[:, imin:imax])), axis=-1)
    wavemax = wave[imin:imax][maxratioloc]

    par = pars[0]

    selection = get_selections(par, sfh=sfh)
        
    fig, axes = pl.subplots(3, 2, figsize=(13, 15))
    for i, (sel, (inds, axnum)) in enumerate(selection.iteritems()):
        ax = axes.flat[axnum]
        for iwave in ishow:
            ax.plot(sratio[inds, iwave], 'o', alpha=0.5, label='$\lambda={}$'.format(wave[iwave]))
        ax.plot(np.exp(maxratio[inds]), 'o', alpha = 0.5, label='max, ${}<\lambda<{}$'.format(wlo, whi))
        ax.plot(mratio[inds], 'o', alpha = 0.5, label='mass')
        ax.set_xlabel('model #')
        ax.set_ylabel('flux(old)/flux(new)')
        ax.set_title(sel)
        #ax.set_ylim(0.7, 1.3)
        ax.axhline(1.0, linestyle=':', color='k')
        if axnum == 0:
            ax.legend(loc=0, prop={'size':10.})
    fig.suptitle('SFH={}'.format(sfh))
    wave = sps.wavelengths
    pl.show()
    return fig, axes, selection, maxratio
    

def get_selections(par, sfh):

    if (sfh == 1) or (sfh == 4):
        simple = ((par['add_neb_emission'] == 0) & (par['dust2'] ==0.0) & (par['fburst'] == 0) &
                  (par['const'] == 0) & (par['sf_start'] < par['sf_trunc'])
                  )

        lesssimple = ((par['add_neb_emission'] == 0) & (par['dust2'] ==0.0) &
                      (par['sf_start'] < par['sf_trunc'])
                     )

        dustysimple = ( (par['add_neb_emission'] == 0) & (par['fburst'] == 0) & (par['const'] == 0) &
                       (par['sf_start'] < par['sf_trunc'])
                     )

        nebsimple = (  (par['dust2'] ==0.0) & (par['fburst'] == 0) & (par['const'] == 0) &
                       (par['sf_start'] < par['sf_trunc'])
                     )

        complicated =  ((par['sf_start'] < par['sf_trunc'])
                        )

        return {'simplest': (simple, 0),
                 'const or fburst': (lesssimple, 1),
                 'dust': (dustysimple, 2),
                 'neb': (nebsimple, 3),
                 'dust, neb, const or fburst': (complicated, 4)
                 }

    if sfh == 5:
        simple = (((par['sf_trunc'] >= par['tage']) | (par['sf_trunc'] == 0)) &
                  (par['add_neb_emission'] == 0))

        lesssimple = (((par['sf_trunc'] < par['tage']) & (par['sf_trunc'] > 0)) &
                      (par['add_neb_emission'] == 0) & (par['sf_slope'] == 0))

        nebsimple = ((par['sf_trunc'] > par['tage']) | (par['sf_trunc'] == 0))

        hislope =  (((par['sf_trunc'] < par['tage']) & (par['sf_trunc'] > 0)) &
                      (par['add_neb_emission'] == 0) & (par['sf_slope'] > 0))

        loslope =  (((par['sf_trunc'] < par['tage']) & (par['sf_trunc'] > 0)) &
                      (par['add_neb_emission'] == 0) & (par['sf_slope'] < 0))

        return {'no truncation': (simple, 0),
                'sf_slope=0': (lesssimple, 1),
                #'neb': (nebsimple, 2),
                'sf_slope = 0.5': (hislope, 2),
                'sf_slope < -0.5': (loslope, 3),
                }

    
if __name__ == "__main__":
    fig, ax, sel, maxratio = main(sfh, root)

    try:
        par = pars[0]
        lesssimple = sel['const or burst']
        bad = lesssimple & (np.exp(maxratio) < 1.5) & (np.exp(maxratio) > 1.08)
        question = lesssimple & (par['tage']-par['sf_start'] > 4) & (par['fburst'] > 0.0)
    except:
        pass
