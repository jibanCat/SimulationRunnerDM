"""
Loading derived summary statistics from nbodykit
from multiple simulations.
"""
from typing import Tuple

import numpy as np

import nbodykit
from nbodykit.lab import ArrayCatalog, FFTPower
import bigfile

def load_nbodykit_power(path: str, scale_factor: int, k_max = None,
        subtract_shotnoise: bool = True, times_kcubic: bool = True) -> Tuple[float, np.ndarray]:
    """
    Parameters:
    ----
    path: path to the PART folder
    scale_factor: for double checking
    """
    bigf = bigfile.File(path)

    header = bigf.open("Header")

    boxsize = header.attrs['BoxSize'][0]
    redshift = 1./header.attrs['Time'][0] - 1
    a = header.attrs['Time'][0]
    assert a == scale_factor

    # npart : number of particle per side
    Ng = header.attrs['TotNumPart'][1] ** (1/3)
    Ng = int(np.rint(Ng))

    pid_ = bigf.open('1/ID')[:] - 1   # so that particle id starts from 0
    pos_ = bigf.open('1/Position')[:]

    f = ArrayCatalog({'Position': pos_ * 0.001})

    k_mean_particle = 2 * np.pi / (boxsize * 0.001) * Ng 
    k_max = 2 * k_mean_particle

    # compute the power spectrum
    mesh = f.to_mesh(resampler='cic', Nmesh=Ng, position='Position', BoxSize=boxsize*0.001)
    rr = FFTPower(mesh,mode='1d', kmax=k_max)
    Pk = rr.power

    # k modes
    ps = Pk['power'].real - Pk.attrs['shotnoise']

    k0, Pk = Pk['k'],ps*Pk['k']*Pk['k']*Pk['k']/(2*np.pi*np.pi)

