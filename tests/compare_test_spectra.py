import pickle, sys
import numpy as np
import matplotlib.pyplot as pl
import fsps
sps = fsps.StellarPopulation(zcontinuous=1)
wave = sps.wavelengths

sfh = 4
branches = ['localvars', 'new_compsp']
spec, parstruct, masses = [], [], []
for branch in branches:
    filename = 'pickles/sfh{}_{}.pkl'.format(sfh, branch)
    with open(filename, 'rb') as f:
        blob = pickle.load(f)
    branch_pars, branch_spec, m, pardicts, default_sfh = blob
    spec.append(branch_spec)
    parstruct.append(branch_pars)
    masses.append(m)
    
missing = [np.all(branch_spec == 0., axis=-1) for branch_spec in spec]
either_missing = missing[0] | missing[1]
both_missing = missing[0] & missing[1]

# These are places where sf_trunc < sf_start, so there should be no SF.
# But new_compsp computes a spectrum!
diff_missing = missing[0] & ~missing[1]
npresent = (~either_missing).sum()

ishow = [100, 500, 700]

wlo, whi = 1e3, 2e4
imin, imax = np.argmin(np.abs(wave - wlo)), np.argmin(np.abs(wave - whi)) 

sdiff = ((spec[0] - spec[1]) / spec[0])
sratio = spec[0]/ spec[1]
mratio = masses[0] / masses[1]
maxratio = np.max(np.abs(np.log(sratio[:, imin:imax])), axis=-1)
#maxratioloc = np.argmax(np.abs(np.log(sratio)), axis=-1)

par = parstruct[0]
simple = ((par['add_neb_emission'] == 0) & (par['dust2'] ==0.0) & (par['fburst'] == 0) &
          (par['const'] == 0) & (~either_missing) & (par['sf_start'] < par['sf_trunc'])
          )

lesssimple = ( (par['add_neb_emission'] == 0) & (par['dust2'] ==0.0) & (~either_missing) &
               (par['sf_start'] < par['sf_trunc'])
             )

dustysimple = ( (par['add_neb_emission'] == 0) & (par['fburst'] == 0) & (par['const'] == 0) & (~either_missing) &
               (par['sf_start'] < par['sf_trunc'])
             )

nebsimple = (  (par['dust2'] ==0.0) & (par['fburst'] == 0) & (par['const'] == 0) & (~either_missing) &
               (par['sf_start'] < par['sf_trunc'])
             )
    
complicated =  ((~either_missing) &
               (par['sf_start'] < par['sf_trunc'])
                )

selection = {'simplest': (simple, 0),
             'const or fburst': (lesssimple, 1),
             'dust': (dustysimple, 2),
             'neb': (nebsimple, 3),
             'dust, neb, const or fburst': (complicated, 4)
             }
    
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

sfh = 4
branches = ['localvars', 'new_compsp']
spec, parstruct = [], []
for branch in branches:
    filename = 'pickles/sfh{}_{}.pkl'.format(sfh, branch)
    with open(filename, 'rb') as f:
        blob = pickle.load(f)
    branch_pars, branch_spec, m, pardicts, default_sfh = blob
    spec.append(branch_spec)
    parstruct.append(branch_pars)
    
missing = [np.all(branch_spec == 0., axis=-1) for branch_spec in spec]
either_missing = missing[0] | missing[1]
both_missing = missing[0] & missing[1]

# These are places where sf_trunc < sf_start, so there should be no SF.
# But new_compsp computes a spectrum!
diff_missing = missing[0] & ~missing[1]
npresent = (~either_missing).sum()

ishow = [100, 500, 700]
pl.show()
bad = lesssimple & (np.exp(maxratio) < 1.5) & (np.exp(maxratio) > 1.08)
question = lesssimple & (par['tage']-par['sf_start'] > 4) & (par['fburst'] > 0.0)
sys.exit()


# Old way of doing things
maxdiff = np.abs(sdiff).max(axis=-1)
wheremaxdiff = np.argmax(np.abs(sdiff), axis=-1)
bad =  (maxdiff > 0.1) & (~either_missing)

#plot(maxdiff(bad))

# these are tage=6.73, sf_start=6.33, sf_trunc=[0, 10].
# localvars returns a spectrum of 1e-70
# 192 models
# Only show up with MILES, time_res_inc=1
superbad = (maxdiff > 2) & (~either_missing)

# 408 models.  These are tau=0.1, tage-sf_start <= 1.0
# new_compsp returns a spectrum of 1e-70, and a mass of nan
badone = (maxdiff == 1.0) & (~either_missing)

# 512 models.  mostly bad in the ionizing flux
bad = (maxdiff != 1.0) & (maxdiff < 2) & (maxdiff > 0.2) & (~either_missing)

# mostly due to nebular lines or ionizing continuum -> young stars
# Some due to NIR flux = tpagb?
#~750 models
# ~
littlebad = (maxdiff < 0.2) & (maxdiff > 0.02) & (~either_missing)

# ~20% less NIR flux in new_compsp for tage=10, sf_start=9
nirbad = (wheremaxdiff == 4841)

# 148 models
# !only 16 with basel, time_resinc=10!  all dust2=1
good = (maxdiff < 0.02) & (~either_missing)
