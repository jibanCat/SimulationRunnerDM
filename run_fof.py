"""
Compute missing FOF tables (if needed)
"""
import os
import numpy as np

from SimulationRunner.multi_haloes import HaloMassFunction

def run_fof(submission_dir: str = "test/", gadget_dir: str = "~/codes/MP-Gadget/",
        mpgadget_param_file: str = "mpgadget.param") -> None:
    """
    Check if there are missing FOF tables. If yes, write a submission file.

    Parameters:
    ----
    submission_dir : the folder includes the MP-Gadget submission files.
    gadget_dir : the absolute path to the MP-Gadget folder.
    mpgadget_param_file : the filename of the MP-Gadget parameter file.
    """
    hmfloader = HaloMassFunction(submission_dir, gadget_dir)

    # if no metal return in the parameter file (old MP-Gadget)
    # add a line to the param file; the returned param_file will
    # be a different file
    mpgadget_param_file = check_metal_return(
        mpgaget_param_path=os.path.join(submission_dir, gadget_dir),
        metal_return=0, # force it to be 0 here since dm-only
    )
    print("[Info] Current param file is {}".format(mpgadget_param_file))

    hmfloader.make_foftables(mpgadget_param_file=mpgadget_param_file)


def check_metal_return(
        mpgadget_param_path: str = "mpgadget.param",
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
    with open(mpgadget_param_path, "r") as f:
        txt = f.readlines()

    if not np.any(["MetalReturnOn" in line for line in txt]):
        txt.append("MetalReturnOn = {}\n".format(metal_return))

        new_file = "metal." + mpgadget_param_path

        with open(new_file, "w") as f:
            f.write("".join(txt))

        return new_file
    
    # if already has MetalReturnOn
    return mpgadget_param_path

