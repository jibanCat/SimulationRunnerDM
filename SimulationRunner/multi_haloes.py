"""
Compute halo mass functions and save into a HDF5 file
"""
from typing import Generator, List

import os
import numpy as np
import re

from .multi_sims import GadgetLoad
from .hmffromfof import HMFFromFOF

class MPISubmit(object):
    
    """
    Read the mpi_submit, the sbatch submission for MP-Gadget
    into a class.
    """

    def __init__(self, filename: str, gadget_dir: str = "~/codes/MP-Gadget/"):
        self.filename = filename
        self:gadget_dir = gadget_dir

        with open(self.filename, "r") as f:
            self.txt = f.readlines()

        
    def get_basic_setups(self):
        """
        Return the lines without `mpirun` - this (ideally) includes
        - SBATCH parameters
        - module loads
        - exports
        """
        # only check the 1st position in a line;
        # this is because sbatch job's name could have mpirun
        exclusion_bag = {"mpirun"}

        basic_txt = []

        for line in self.txt:
            words = line.split()

            if np.any([words[0] == exclusion_word for exclusion_word in exclusion_bag]):
                continue
    
            basic_txt.append(line)

        return basic_txt

    def make_simulation_foftable(self, snapshots: List[int], mpi_submit_file: str = "mpi_submit_foftables") -> None:
        """
        Generate a submission file for making the fof table (PART -> PIG)
        """
        basic_txt = self.get_basic_setups()

        # snapshots could be more than one
        for snapshot in snapshots:
            basic_txt.append(self.mpirun_foftable(snapshot))
        

        with open(mpi_submit_file, "r") as f:
            f.write("".join(basic_txt))


    @property
    def mpirun_foftable(self, snapshot: int) -> str:
        """
        The line to make foftable from PART snapshot       
        """
        return "mpirun --map-by core {} mpgadget.param 3 {}\n".format(self.gadget_dir, snapshot)

class HaloMassFunction(GadgetLoad):

    """
    Compute the halo mass functions from foftable using MP-Gadget's
    tool, HMFFromFOF (TODO: binning scheme might need to adjust,
    https://jakevdp.github.io/blog/2012/09/12/dynamic-programming-in-python/)

    There was a bug before that FOF tables were not generated from
    the PART during the run time. Generate a script to get FOF tables
    out of PART.

    Example:
    ./MP-Gadget paramfile.gadget 3 $snapnum
    """

    def __init__(self, submission_dir: str = "test/") -> None:
        super(HaloMassFunction, self).__init__(submission_dir)

        # acquire the current PART files
        self._parts = sorted(self.get_PART_array())
        # acquire the current PIG files
        self._pigs = sorted(self.get_PIG_array())

        # read mpi_submit into a class
        # keep all SBATCH parameters the same;
        # keep all module load the same
        # make the mpirun free parameter

        # check # of FOF tables == # of snapshots
        self._pigs_to_run = set(self._parts) - set(self._pigs)
        if len(self._pigs_to_run) > 0:
            print("[Warning] Some snapshots lack of FOF tables, recommend running ...")

    def make_foftables(self, snapshots: List[int]):
        """
        Make FOF tables from PART
        """




    def get_PART_array(self) -> Generator:
        """
        Get indexes for PART snapshot folders
        """
        for f in self._outputfiles:
            out = re.findall(r"PART_([0-9]*)", f)

            # if matched
            if len(out) == 1:
                yield int(out[0])
    
    def get_PIG_array(self) -> Generator:
        """
        Get indexes for PIG snapshot folders
        """
        for f in self._outputfiles:
            out = re.findall(r"PIG_([0-9]*)", f)

            # if matched
            if len(out) == 1:
                yield int(out[0])
    




