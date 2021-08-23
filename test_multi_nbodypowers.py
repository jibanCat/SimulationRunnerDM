"""
Test file for creating hdf5 catalogue for powerspecs from multiple simulations,
using the powerspecs computed from nbodykit
"""
from typing import List, Optional
import os
import numpy as np
import h5py
from SimulationRunner.multi_sims import fn_outdir
from SimulationRunner.multi_nbodykit import HDF5Holder, MultiNbodyKitPowerSpec, NbodyKitPowerSpec, HDF5Holder

from test_multi_powerspecs import test_fileload

def test_customRes(
    base_dir: str, res: int, box: int, n_simulations: int, Latin_json: str,
    srgan: bool = False, z0 : float = 0.0, Ng: int = 512, kmax: float = 16.10,
    srgan_path: str = "super-resl/output/PART_008/powerspec_shotnoise.txt.npy",
):
    """
    Create a custom Latin Hyper Cube multi-spec file.
    """
    # This is just creating a list of filenames
    test_dir = lambda i: os.path.join(base_dir, fn_outdir(i, res, box))

    all_submission_dirs = [test_dir(i) for i in range(n_simulations)]

    # test lowRes hdf5
    test_create_hdf5(
        all_submission_dirs, Latin_json=os.path.join(base_dir, Latin_json),
        srgan=srgan, z0=z0, Ng=Ng, kmax=kmax, srgan_path=srgan_path
    )


def test_create_hdf5(
    all_submission_dirs: List[str],
    Latin_json: str,
    srgan: bool, z0 : float, Ng: int, kmax: float,
    srgan_path: str,
    selected_ind: Optional[np.ndarray] = None,
    ) -> None:
    multinps = MultiNbodyKitPowerSpec(
        all_submission_dirs, Latin_json, selected_ind=selected_ind,
        srgan=srgan, z0=z0, Ng=Ng, kmax=kmax, srgan_path=srgan_path)

    # test you have simulations in the dirs
    test_fileload(multinps.all_submission_dirs[0])
    test_fileload(multinps.all_submission_dirs[-1])

    # sample size should be the same
    assert len(multinps.Latin_dict["omega0"]) == len(multinps.Latin_dict["omegab"])
    assert len(multinps.Latin_dict["hubble"]) == len(multinps.Latin_dict["omegab"])

    multinps.create_hdf5("test_dmonly.hdf5")

    assert os.path.exists("test_dmonly.hdf5")

    with h5py.File("test_dmonly.hdf5", "r") as test_hdf5:
        # test the first simulation
        assert "simulation_0" in test_hdf5
        assert "simulation_0/powerspecs" in test_hdf5
        assert "simulation_0/k0" in test_hdf5
        if srgan:
            assert "simulation_0/powerspecs_srgan" in test_hdf5

        # we condition on z0=0
        redend = test_hdf5["simulation_0"].attrs["redend"]
        scale = 1 / (redend + 1)

        sim = test_hdf5["simulation_0"]
        assert np.abs(scale - sim["scale_factors"][-1]) < 1e-6

        # make sure the z0 is the same
        assert (1 / (z0 + 1)) == sim["scale_factors"][0]

        # test the ks SRGAN: should be the same length as powerspecs
        if srgan:
            assert np.all(sim["k0_sr"][()] == sim["k0"][()])

        # test the sampling size is correct
        length = len(test_hdf5["hubble"]) - 1

        assert "simulation_{}".format(length) in test_hdf5
        assert "simulation_{}/powerspecs".format(length) in test_hdf5

        # test the ordering of Latin HyperCube is the same as simulation_{}
        hubble1 = test_hdf5["hubble"][0]
        hubble11 = test_hdf5["simulation_0"].attrs["hubble"]

        assert np.abs(hubble1 - hubble11) < 1e-6

        hubble_end = test_hdf5["hubble"][-1]
        hubble_endend = test_hdf5["simulation_{}".format(length)].attrs["hubble"]

        assert np.abs(hubble_end - hubble_endend) < 1e-6

    # output the txt files
    f = HDF5Holder("test_dmonly.hdf5")
    f.to_txt(srgan_output=srgan)
