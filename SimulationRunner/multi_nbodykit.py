"""
Loading derived summary statistics from nbodykit
from multiple simulations.
"""
from SimulationRunner.multi_sims import PowerSpec
from typing import Tuple, List

import os

import numpy as np

import nbodykit
from nbodykit.lab import ArrayCatalog, FFTPower
import bigfile

def load_nbodykit_power(path: str, scale_factor: int, k_max = None,
        subtract_shotnoise: bool = True, times_kcubic: bool = False, compensated=True,
        Ng=None,) -> Tuple[float, np.ndarray]:
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
    if Ng == None:
        Ng = header.attrs['TotNumPart'][1] ** (1/3)
        Ng = int(np.rint(Ng))

    pid_ = bigf.open('1/ID')[:] - 1   # so that particle id starts from 0
    pos_ = bigf.open('1/Position')[:]

    f = ArrayCatalog({'Position': pos_ * 0.001})

    # compute until 2 times mean particle spacing if k_max not given
    if k_max == None:
        k_mean_particle = 2 * np.pi / (boxsize * 0.001) * Ng 
        k_max = 2 * k_mean_particle

    # compute the power spectrum
    mesh = f.to_mesh(resampler='cic', Nmesh=Ng, position='Position', BoxSize=boxsize*0.001, compensated=compensated)
    rr = FFTPower(mesh,mode='1d', kmax=k_max)
    Pk = rr.power

    k0 = Pk["k"]

    # original power
    ps = Pk['power'].real

    if subtract_shotnoise:
        ps = Pk['power'].real - Pk.attrs['shotnoise']
    
    if times_kcubic:
        ps = Pk['power'].real - Pk.attrs['shotnoise']
        ps = ps*Pk['k']*Pk['k']*Pk['k']/(2*np.pi*np.pi)
    
    return k0, ps

class NbodyKitPowerSpec(PowerSpec):

    """
    Loading power spectra from Nbodykit.

    Loading all snaphots for MP-Gadget outputs.
    Loading one snapshot for SR output.

    Note: the training set will be directly outputed from MultiNbodyPowerSpec
    """

    def __init__(self, 
        submission_dir: str = "test/", srgan: bool = False, z0 : float = 0.0, Ng: int = 512,
        srgan_path: str = "super-resl/output/PART_008/powerspec_shotnoise.txt.npy") -> None:
        super(NbodyKitPowerSpec, self).__init__(submission_dir)

        # read into arrays
        # Matter power specs from simulations
        k0, ps = self.read_powerspec(z0=z0, Ng=Ng)

        self._scale_factors = 1 / (1 + z0)

        self._k0 = k0
        self._powerspecs = ps

        # Matter power specs from CAMB linear theory code
        redshifts, out = self.read_camblinear(self.camb_files)

        self._camb_redshifts = redshifts
        self._camb_matters = out

        # Matter power specs from SRGAN (conditioned on one redshift)
        if srgan:
            self.srgan_path = srgan_path
            self._k0_sr, self._ps_sr = self.read_srgan_powerspec(srgan_path)

    @property
    def powerspecs(self) -> np.ndarray:
        """
        P(k) from PART/ folder. The same length as k0.
        """
        return self._powerspecs

    @property
    def k0(self) -> np.ndarray:
        """
        k from PART/ folder. The same length as powerspecs
        """
        return self._k0

    @property
    def powerspecs_srgan(self) -> np.ndarray:
        """
        P(k) from SRGAN.
        """
        return self._ps_sr
    
    @property
    def k0_sr(self) -> np.ndarray:
        """
        k for SRGAN power spectrum.
        """
        return self._k0_sr

    def read_powerspec(self, z0: float, Ng: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Read power spectrum from a PART/ folder

        Parameters:
        ----
        z0 : the redshift of the power spectrum you want to load.
        Ng : the number of particle per side you want to compute for the power spectrum.
        """
        tol = 1e-4 # tolerance
        # the scale factor you condition on
        scale_factor = 1 / (1 + z0)

        # | No. of snapshot | scale factor |
        # self._snapshots
        if not np.any((self._snapshots[:, 1] - scale_factor) < tol):
            print("Pick a scale factor in the list:", self._snapshots[:, 1])
            assert np.any((self._snapshots[:, 1] - scale_factor) < tol)

        # -1: if there are multiple matches, pick the final one (just personal preference)
        ii = np.where((self._snapshots[:, 1] - scale_factor) < tol)[0][-1]

        # PART_{number}
        print("Found a={} is PART_{}.".format(scale_factor, self._snapshots[:, 0][ii]))
        number = int(self._snapshots[:, 0][ii])

        powerspec_path = os.path.join(
            self.submission_dir, "output", "PART_{:03d}".format(number)
        )

        # the maximum k is controlled by Ng
        k0, ps = load_nbodykit_power(
            powerspec_path, scale_factor=scale_factor, k_max=None, subtract_shotnoise=False, times_kcubic=False, compensated=True, Ng=Ng)

        return k0, ps


    def read_srgan_powerspec(self, srgan_path: str) -> Tuple[np.ndarray, np.ndarray]:
        """
        Read the SRGAN power spectrum directly from the file.

        SRGAN is conditioned on a single redshift, so be aware which redshift you applied before.
        """
        # load the SR
        # exception happens when I wrongly saved the power spec as a binary file
        try:
            k0_sr, ps_sr = np.loadtxt(srgan_path)
        except UnicodeDecodeError as e:
            k0_sr, ps_sr = np.load(srgan_path, allow_pickle=True)
        
        return k0_sr, ps_sr


class MultiNbodyKitPowerSpec:
    """
    Output a HDF5 file.

    An additional function to output the training power spectra directly,
    includes:
        1. interpolate the LowRes (if necessary).
        2. substract surrogate mean.
        3. condition on z = 0.

    """
    def __init__(self) -> None:
        pass