import numpy as np
import matplotlib.pyplot as pl


with open('../data/FILTER_LIST', 'r') as f:
    filts = f.readlines()
fnames = [fi.split()[1].lower() for fi in filts]

def read_mags(branch):
    fn = '{}_CSP_tabsfh_test.out.mags'.format(branch)
    with open(fn, 'r') as f:
        lines = f.readlines()
    hdr = lines[7][1:].split()
    hdr = hdr[:4] + fnames

    mags = np.genfromtxt(fn, skiprows=8, dtype=np.dtype(zip(hdr, len(hdr) * [np.float])))
    return mags

def read_sfh():
    fn = '../data/sfh.dat'
    age, sfr, z = np.genfromtxt(fn, unpack=True, skiprows=0)
    return age, sfr, z
    
if __name__ == "__main__":
    oldmags = read_mags('localvars')
    newmags = read_mags('newcompsp')
    ages, sfrs, z = read_sfh()

    
    fig, axes = pl.subplots(2, 1, sharex=True)
    ax = axes[0]
    ax.plot(oldmags['logage'], oldmags['v'], '-o', label='old')
    ax.plot(newmags['logage'], newmags['v'], '-o', label='new')
    ax.set_ylim(-16, -25)

    ax = axes[1]
    ax.plot(oldmags['logage'], oldmags['logSFR'], '-o', label='old')
    ax.plot(newmags['logage'], newmags['logSFR'], '-o', label='new')
    ax.plot(np.log10(ages) + 9, np.log10(sfrs), '-o', label='input')
    ax.set_ylim(-1, 2)

    [a.set_xlim(8.5, 10.3) for a in axes]
    [a.legend(loc=0) for a in axes]
    fig.show()
