'''
Test file for creating hdf5 catalogue from multiple simulations
'''
from typing import List
import os
import numpy as np
import h5py
from SimulationRunner.multi_sims import powerspec_fn, MultiPowerSpec

def test_fileload(submission_dir : str) -> None:
    # test if you have SimulationICs.json
    simics = os.path.join(submission_dir, "SimulationICs.json")

    assert os.path.exists(simics)

    # test if you run until a = 1
    powerspec = os.path.join(submission_dir, "output", powerspec_fn(1))

    assert os.path.exists(powerspec)

def test_create_hdf5(all_submission_dirs: List[str], Latin_json: str):
    multips = MultiPowerSpec(all_submission_dirs, Latin_json)

    # test you have simulations in the dirs
    test_fileload(multips.all_submission_dirs[0])
    test_fileload(multips.all_submission_dirs[-1])

    # sample size should be the same
    assert len(multips.Latin_dict['omega0']) == len(multips.Latin_dict['omegab'])
    assert len(multips.Latin_dict['hubble']) == len(multips.Latin_dict['omegab'])

    multips.create_hdf5("test_dmonly.hdf5")

    assert os.path.exists("test_dmonly.hdf5")

    with h5py.File("test_dmonly.hdf5", "r") as test_hdf5:
        # test the first simulation
        assert "simulation_1" in test_hdf5
        assert "simulation_1/powerspecs" in test_hdf5
        
        redend = test_hdf5["simulation_1"].attrs["redend"]
        scale  = 1 / (redend + 1)
        
        sim = test_hdf5["simulation_1"]
        assert np.abs(scale - sim["scale_factors"][-1]) < 1e-6

        # test the sampling size is correct
        length = len(test_hdf5["hubble"]) - 1

        assert "simulation_{}".format(length) in test_hdf5
        assert "simulation_{}/powerspecs".format(length) in test_hdf5

        # test the ordering of Latin HyperCube is the same as simulation_{}
        hubble1  = test_hdf5["hubble"][0]
        hubble11 = test_hdf5["simulation_1"].attrs["hubble"]

        assert np.abs(hubble1 - hubble11) < 1e-6

        hubble_end    = test_hdf5["hubble"][-1]
        hubble_endend = test_hdf5["simulation_{}".format(length)].attrs["hubble"]

        assert np.abs(hubble_end - hubble_endend) < 1e-6