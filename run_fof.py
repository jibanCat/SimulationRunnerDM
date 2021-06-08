"""
Compute missing FOF tables (if needed)
"""
import os
import numpy as np

from SimulationRunner.multi_haloes import HaloMassFunction

def run_fof(submission_dir: str = "test/", gadget_dir: str = "~/codes/MP-Gadget/") -> None:
    """
    Check if there are missing FOF tables. If yes, write a submission file.
    """
    hmfloader = HaloMassFunction(submission_dir, gadget_dir)



def check_metal_return(
        mpgaget_param_path: str = "mpgadget.param",
        metal_return: int = 0
    ) -> str:
    """
    Check if the parameter file has MetalReturnOn.

    Parameters:
    ----
    metal_return: MetalReturnOn = 0 or 1, Enable the return of metals from star particles to the gas

    Return:
    ---
    submission_file_name
    """
    with open(mpgaget_param_path, "r") as f:
        txt = f.readlines()

    if not np.any(["MetalReturnOn" in line for line in txt]):
        txt.append("MetalReturnOn = {}\n".format(metal_return))

        new_file = "metal." + mpgaget_param_path

        with open(new_file, "w") as f:
            f.write("".join(txt))

        return new_file
    
    # if already has MetalReturnOn
    return mpgaget_param_path

