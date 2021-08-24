"""
Loading derived summary statistics from nbodykit
from multiple simulations.
"""
from SimulationRunner.multi_sims import MultiPowerSpec, PowerSpec
from typing import Tuple, List, Optional, Generator

import os

import numpy as np
import h5py

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
        submission_dir: str = "test/", srgan: bool = False, z0 : float = 0.0, Ng: int = 512, kmax: float =16.10,
        srgan_path: str = "super-resl/output/PART_008/powerspec_shotnoise.txt.npy") -> None:
        super(PowerSpec, self).__init__(submission_dir)

        # read into arrays
        # Matter power specs from simulations
        k0, ps = self.read_powerspec(z0=z0, Ng=Ng, kmax=kmax)

        self._scale_factors = np.array([1 / (1 + z0)])

        self._k0 = k0
        self._powerspecs = ps

        # Matter power specs from CAMB linear theory code
        redshifts, out = self.read_camblinear(self.camb_files)

        self._camb_redshifts = redshifts
        self._camb_matters = out

        # Matter power specs from SRGAN (conditioned on one redshift)
        if srgan:
            self.srgan_path = os.path.join(submission_dir, srgan_path)
            self._k0_sr, self._ps_sr = self.read_srgan_powerspec(self.srgan_path, kmax=kmax)

            # [TODO] you need write a interpolate function if they don't match
            assert np.all(self._k0_sr == self._k0)

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

    def read_powerspec(self, z0: float, Ng: int, kmax: float) -> Tuple[np.ndarray, np.ndarray]:
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

        # filter out NaN values
        ind = ~np.isnan(k0)
        assert np.all(ind == ~np.isnan(ps))

        # set the kmax
        ind = ind & (k0 <= kmax)
        # remove k = 0 since no power there
        ind = ind & (k0 != 0.0)

        k0 = k0[ind]
        ps = ps[ind]

        return k0, ps


    def read_srgan_powerspec(self, srgan_path: str, kmax: float, ) -> Tuple[np.ndarray, np.ndarray]:
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

        # filter out NaN values
        ind = ~np.isnan(k0_sr)
        assert np.all(ind == ~np.isnan(ps_sr))
        # set the kmax
        ind = ind & (k0_sr <= kmax)
        # remove k = 0 since no power there
        ind = ind & (k0_sr != 0.0)

        k0_sr = k0_sr[ind]
        ps_sr = ps_sr[ind]

        return k0_sr, ps_sr


class MultiNbodyKitPowerSpec(MultiPowerSpec):
    """
    Output a HDF5 file.

    An additional function to output the training power spectra directly,
    includes:
        1. interpolate the LowRes (if necessary).
        2. substract surrogate mean.
        3. condition on z = 0.

    """
    def __init__(self, all_submission_dirs: List[str], Latin_json: str, selected_ind: Optional[np.ndarray],
        srgan: bool = False, z0 : float = 0.0, Ng: int = 512, kmax: float =16.10,
        srgan_path: str = "super-resl/output/PART_008/powerspec_shotnoise.txt.npy") -> None:
        super().__init__(all_submission_dirs, Latin_json=Latin_json, selected_ind=selected_ind)

        # assign attrs for loading Nbodykit power specs
        self.srgan = srgan
        self.z0 = z0
        self.Ng = Ng
        self.srgan_path = srgan_path
        self.kmax = kmax


    def create_hdf5(self, hdf5_name: str = "MutliPowerSpecs.hdf5") -> None:
        """
        - Create a HDF5 file for powerspecs from multiple simulations.
        - Each simulation stored in subgroup, includeing powerspecs and
        camb linear power specs.
        - Each subgroup has their own simulation parameters extracted from
        SimulationICs.json to reproduce this simulation.
        - Parameters from Latin HyperCube sampling stored in upper group level,
        the order of the sampling is the same as the order of simulations.

        TODO: add a method to append new simulations to a created hdf5.
        """
        # open a hdf5 file to store simulations
        with h5py.File(hdf5_name, "w") as f:
            # store the sampling from Latin Hyper cube dict into datasets:
            # since the sampling size should be arbitrary, we should use
            # datasets instead of attrs to stores these sampling arrays
            for key, val in self.Latin_dict.items():
                f.create_dataset(key, data=val)

            # using generator to iterate through simulations,
            # PowerSpec stores big arrays so we don't want to load
            # everything to memory
            for i, ps in enumerate(self.load_PowerSpecs(self.all_submission_dirs, srgan=self.srgan, z0=self.z0, Ng=self.Ng, kmax=self.kmax, srgan_path=self.srgan_path)):
                sim = f.create_group("simulation_{}".format(i))

                # store arrays to sim subgroup
                # Simulation Power spectra:
                sim.create_dataset("scale_factors", data=np.array(ps.scale_factors))
                sim.create_dataset("powerspecs", data=ps.powerspecs)
                sim.create_dataset("k0", data=ps.k0)
                # SRGAN power spectra:
                if self.srgan:
                    sim.create_dataset("powerspecs_srgan", data=ps.powerspecs_srgan)
                    sim.create_dataset("k0_sr", data=ps.k0_sr)

                sim.create_dataset("camb_redshifts", data=np.array(ps.camb_redshifts))
                sim.create_dataset("camb_matters", data=ps.camb_matters)

                # stores param json to metadata attrs
                for key, val in ps.param_dict.items():
                    sim.attrs[key] = val

    @staticmethod
    def load_PowerSpecs(all_submission_dirs: List[str], srgan: bool, z0 : float, Ng: int, kmax: float,
            srgan_path: str) -> Generator:
        """
        Iteratively load the PowerSpec class
        """
        for submission_dir in all_submission_dirs:
            yield NbodyKitPowerSpec(submission_dir, srgan=srgan, z0=z0, Ng=Ng, kmax=kmax, srgan_path=srgan_path)

class HDF5Holder(h5py.File):
    """
    Hold the h5 file generated by SimulationRunner.multi_sims.MultiPowerSpes,
    and add some class manipulation methods for it.
    """

    def __init__(self, name: str, mode: str = "r", saved_filename: str = "test.h5"):
        super().__init__(name, mode=mode)

        # a temp file saved in the root
        self.saved_filename = saved_filename
        self._mode = mode

    def interpolate(self, ks: np.ndarray):
        """
        interpolate the log P(k) based on a given ks 
        """
        pass

    def to_txt(self, srgan_output: bool, input_filename: str = "input.txt", output_filename: str = "outupt.txt") -> None:
        """
        output the data to .txt file with:
        ----
        input.txt (number of simulations, number of parameters) : cosmological parameters
        output.txt (number of simulations, number of k modes) : power spectra
        
        The output files will be loaded by matter_multi_fidelity_emu.data_loader.PowerSpecs

        Parameters:
        ----
        srgan_output : use SRGAN's output as output power spectra
        """
        # prepare input parameters
        X = np.array([self[name][()] for name in self["parameter_names"]]).T

        num_simulations, _ = X.shape

        # prepare output power spectra
        if srgan_output: 
            Y = np.stack([self["simulation_{}".format(i)]["powerspecs_srgan"][()] for i in range(num_simulations)])
            k0 = self["simulation_0"]["k0_sr"][()]
            assert np.all(k0 == self["simulation_0"]["k0"][()])
        else:
            Y = np.stack([self["simulation_{}".format(i)]["powerspecs"][()] for i in range(num_simulations)])
            k0 = self["simulation_0"]["k0"][()]

        np.savetxt(input_filename, X)
        np.savetxt(output_filename, np.log10(Y))

        np.savetxt("input_limits.txt", self["bounds"])
        np.savetxt("kf.txt", np.log10(k0))

    def __add__(self, other):
        """
        combine two HDF5 files for multi power spectra
        """
        if os.path.exists(self.saved_filename):
            os.remove(self.saved_filename)

        # make sure the final size
        parameter_names = self["parameter_names"][()]
        assert np.all(parameter_names == other["parameter_names"][()])
        self_size = self[parameter_names[0]].shape[0]
        other_size = other[parameter_names[0]].shape[0]
        assert other_size == other[parameter_names[-1]].shape[0]

        # it is safer to create another file than modifying the original file
        with h5py.File(self.saved_filename, "w") as new:
            new.create_dataset("parameter_names", data=parameter_names)
            assert np.all(self["bounds"][()] == other["bounds"][()])
            new.create_dataset("bounds", data=self["bounds"][()])

            # append Latin Hypercube parameters; not Latin hypercube anymore
            for param in parameter_names:
                val = np.append(self[param][()], other[param][()])
                assert len(self[param][()].shape) == 1
                assert val.shape[0] == (self_size + other_size)
                new.create_dataset(param, data=val)

            # [npart and boxsize checks]
            assert (
                other["simulation_0"].attrs["npart"]
                == self["simulation_0"].attrs["npart"]
            )
            assert (
                other["simulation_0"].attrs["box"] == self["simulation_0"].attrs["box"]
            )

            # append simultions subgroups
            for i in range(self_size):
                sim = new.create_group("simulation_{}".format(i))

                for key in self["simulation_{}".format(i)].keys():
                    sim.create_dataset(
                        key, data=self["simulation_{}".format(i)][key][()]
                    )

                # append attributes of the subgroup
                for key in self["simulation_{}".format(i)].attrs.keys():
                    sim.attrs[key] = self["simulation_{}".format(i)].attrs[key]

            for i in range(other_size):
                sim = new.create_group("simulation_{}".format(i + self_size))

                for key in other["simulation_{}".format(i)].keys():
                    sim.create_dataset(
                        key, data=other["simulation_{}".format(i)][key][()]
                    )

                # append attributes of the subgroup
                for key in other["simulation_{}".format(i)].attrs.keys():
                    sim.attrs[key] = other["simulation_{}".format(i)].attrs[key]

            assert "simulation_{}".format(self_size + other_size - 1) in new.keys()

        return HDF5Holder(self.saved_filename, mode=self._mode)

