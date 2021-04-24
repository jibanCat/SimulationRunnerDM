# Python code to automate the generation of Gadget3 simulation config files

**Note: this repo now only supports running dark-matter only simulations**.

The base class is Simulation, which creates the config files for a single simulation.

Machine-specific data is implemented with a function which dynamically subclasses the base class.

## How to generate submission files for BIOCluster

### Generate submission files for a batch of simulations

In this section, we will outline how to use the `SimulationICs` class to generate submission files for a batch of simulations for [MP-Gadget](https://github.com/MP-Gadget/MP-Gadget).
We will use `run_dmonly.py` to produce the submission files.

First, make sure you have compiled the MP-Gadget software on the biocluster.
Read the instructions in the page of [MP-Gadget](https://github.com/MP-Gadget/MP-Gadget) to know how to compile the software.

My default path to MP-Gadget is `~/bigdata/codes/MP-Gadget/`.
And also make sure you have your python environment loaded.
There are some required python packages needed to be installed:

- Python 3+
- `configobj`
- `numpy`
- `matplotlib`
- `scipy`
- [`nbodykit`](https://github.com/bccp/nbodykit)
- [`classylss`](https://github.com/nickhand/classylss)

Next, we need a file containing a list of input cosmological parameters for the MP-Gadget simulations.
Here we only support using 5 parameters, `["omega0", "omegab", "hubble", "scalar_amp", "ns"]`

I provide an example file here in the `input_latincube/` folder.
This file has 10 samples from a 5-D Latin hypercube.

To generate the submission files, run:

```bash
cd input_latincube
python ../run_dmonly.py --json_file=matterLatin_high.json --gadget_dir=<path/to/MP-Gadget/> --box=256 --npart=128
```

The `--box` flag indicates your preferred boxsize for the simulation. Here we use 256 Mpc/h.
And the `--npart` flag is the number of dark matter particles per box side you want to put into the simulation box. Here we perform a small size of simulation with 128^3 particles.

After you run the `run_dmonly.py`, if there is no error, you will find several folders in your current directory.
They are named `test-<npart>-<box>-dmonly_<ordering in the json file>/`.
In the above example,
the first folder will be named `test-128-256-dmonly_0000/`.

You will have the following files in the folder:
```
camb_linear   _class_params.ini  ICS             mpi_submit        Options.mk  SimulationICs.json
cambpower.py  _genic_params.ini  mpgadget.param  mpi_submit_genic  output
```
`mpi_submit` and `mpi_submit_genic` are the submission files.
`output` is the folder for the simulation output files.
`SimulationICs.json` is a json file containing all necessary information to reproduce the simulation.

Suppose we want to submit the first simulation.
Do:
```bash
cd test-128-256-dmonly_0000
sbatch mpi_submit_genic     # generate initial conditions for MP-Gadget
sbatch mpi_submit
```

After the simulation finished, you will find the output files in the `output/` folder.
Further information about the output files please check the [MP-Gadget manual](https://www.overleaf.com/read/kzksrgnzhtnh).

For example, `output/powerspectrum-<scale factor>.txt` contains the data for the matter power spectrum.
The scale factor is from 0 ~ 1, where 1 means the current Universe.

### Generate submission files for a single simulation

If you only want to perform a single simulation, do this in Python:
```python
from SimulationRunner import simulationics, clusters

# set the desired input parameters for your simulation
box = 256   # Mpc/h
npart = 128 # 128^3
hubble = 0.7 # Hubble parameter, h, which is H0 / (100 km/s/Mpc)
omega0 = 0.288 # Total matter density at z=0 (includes massive neutrinos and baryons)
omegab = 0.0472 # baryon density.
scalar_amp = 2.427e-9 # A_s at k = 0.05, comparable to the Planck value.
ns = 0.97 # Scalar spectral index
nproc = 256 # Total number of processors
cores = 32 # Number of cores per node

# set the paths
outdir = "test-128-256-dmonly-0000" # Output folder name
gadget_dir = "~/codes/MP-Gadget/" # Your path to MP-Gadget folder
python = "python" # Your path to python binary

# here you run on UCR biocluster
# if you don't want to send emails to my email address, change the email in BIOClass class :)
cluster_class = clusters.BIOClass

Sim = simulationics.SimulationICs(
    redshift = 99, redend = 0,
    box    = box, npart = npart,
    hubble = hubble, omega0 = omega0,
    omegab = omegab, scalar_amp = scalar_amp, ns = ns,
    nproc  = nproc, cores=cores,
    outdir = outdir,
    gadget_dir = gadget_dir,
    python = python,
    cluster_class = cluster_class)

Sim.make_simulation(pkaccuracy=0.07)
assert os.path.exists(outdir)
```

If no errors show up, then you will find the submission files in the `test-128-256-dmonly-0000/` folder.
To run the simulation, do this in bash:
```bash
cd test-128-256-dmonly_0000
sbatch mpi_submit_genic     # generate initial conditions for MP-Gadget
sbatch mpi_submit
```


## How to generate Latin hypercube JSON file

You will need python packages:
- `emukit`
- `pyDOE`

The latin hypercube function is from sbird's lyaemu, which uses rejection-sampled latin hypercubes.

The test file, `test_design.py`, can generate a hypercube with 10 samples with the 5-D LCDM parameters, (h, Omega_0, Omega_b, A_s, n_s).

In case you would like to generate different number of samples (e.g., 50 samples) but with the same bounded uniform prior:
```python
import numpy as np
from latin_design.matter_power_design import MatterDesign

# set the bounds for uniform priors
matter = MatterDesign(
    omegab_bounds=(0.0452, 0.0492), 
    omega0_bounds=(0.268, 0.308), 
    hubble_bounds=(0.65, 0.75), 
    scalar_amp_bounds=(1.5*1e-9, 2.8*1e-9), 
    ns_bounds=(0.9, 0.99))

# number of samples
point_count = 50
matter.save_json(point_count, "matterLatin_50.json")
```
