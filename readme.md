
# Medras

This is a python implementation of the analytic MEDRAS (Mechanistic DNA Repair and Survival) model. This model simulates the response of cells to DNA damage following ionising radiation exposure, including DNA repair kinetics, accuracy of repair, and resulting cell fate. 

Descriptions of the model and the underlying analytic analysis have been published in (McMahon 2016) and (McMahon 2017). This code facilitates the simulation of radiation responses for cells based on a number of phenotypic characteristics, without the requirement for explicit measurement of radiation sensitivity parameters. In particular, the involved parameters are the repair genome size, number of chromosomes, repair capacity, activity of G1 arrest, and cell cycle phase.

This model simulates a range of endpoints, including DNA repair, mutation and chromosome aberration formation, and cell death. These can be queried either as terminal endpoints, or at various times following irradiation. Irradiations can be described as either single instantaneous doses, or protracted exposures potentially with inter-fraction gaps (Note: Cell cycle progression during irradiation is currently not considered).

The key features of the model are briefly described below.

## Using Medras

The core Medras object is a `singleCell` object, which defines the features of a single cell, and provides methods to simulate different types of exposures, defined in the `medrascell.py` file. By default, this simulates a normal human fibroblast in G1, whose survival response can be simulated as follows:

```py
# Add medras to path as preferred
from medras import medrascell

# Generate a simple dose response curve
testCell = medrascell.singleCell() # Initialise cell object
for d in [0,1,2,3,4,5,6]:
	exposure = {'dose':d}
	print(d,'\t',testCell.survival(exposure)) # Print survival for each dose
```

The first command generates a `singleCell` object, populated with default parameters. Radiation responses can then be simulated to different exposures of ionising radiation, the simplest of which can be defined simply in terms of the dose. This is then assumed to be an X-ray irradiation delivered in negligible time. (See below for further details on exposure definition.) This code will thus generate a dose response curve from 0 to 6 Gy for these cells, printed to the terminal.

### Defining cells

Cell characteristics are defined in terms of a dictionary object containing all of the relevant information for a simplified "radiation response phenotype". The default settings for these are:

```py
defaultCell = {'dna':6100.0,'chromosomes':46,'repair':0,'G1Arrest':1,'phase':0,'gene':0}
```

In order, these are:

- DNA content: in MBP, for the full genetic content of the cell (e.g. 6100 MBP for human diploid cells).
- Chromosome number (total);
- Repair defect: stored as a bit mask, 0 for normal cells, 1 for Non-homologous end joining-defective, 2 for homologous recombination-defective, and 3 for double-defect cells;
- G1Arrest status: Activity of the G1 arrest checkpoint. This can be any value between fully active (1) or fully defective (0). This defective checkpoint is heavily associated with p53 activity, with many p53 mutant cells losing this checkpoint;
- Phase: The cell cycle phase of the cell, defined as an integer from 0-3 (from G1 to M, respectively). -1 can be used to simulate an asynchronous population for survival endpoints;
- Gene size: in MBP, this defines the size of a gene of interest for simulating mutation rates, if desired.

A new cell line can be specified varying any or all of these parameters, with the remainder initialised from the defaults. Thus, a normal fibroblast in G2 may be initialised as `singleCell({'phase':2})`, while code to simulate an NHEJ-defective Chinese Hamster Ovary variant in G1 would be:

```py
from medras import medrascell
# Generate a dose response curve for NHEJ-defective hamster
cellParams = {'dna':5400.0,'chromosomes':22,'repair':1,'G1Arrest':0,'phase':0}
testCell = medrascell.singleCell(cellParams) # Initialise alternative cell
for d in [0,1,2,3,4,5,6]:
	exposure = {'dose':d}
	print(d,'\t',testCell.survival(exposure)) # Print survival for each dose
```

In principle any combination of these features should be valid, but it's not possible to guarantee these will be meaningful far from 'typical' mammalian cell characteristics.

### Defining exposures

Much like cells, radiation exposures are defined using dictionaries. The defined features can include:

```py
exposure = {'dose':2, 'time':0, 'LET':5, 'particle':6}
```

In order, these are:

- dose: Delivered to the cells. This can be defined either as a single instantaneous dose in Gy, as here, or a protracted exposure (see below).
- time: This is the time delay applied to the experiment. For measurements of DNA damage, mutations, or aberrations, this is taken as the time at which the measurement is carried out. For survival, this is taken as the plating delay before cells are allowed to proliferate (set to 0 for cells proliferating at the time of irradiation). If a value of -1 is used, the final state of the cell after all repair has completed will be returned. 
- LET: The LET of the incident radiation, in keV/μm. 
- particle: The atomic number of the incident radiation.

Different radiation qualities are handled by combining pre-computed particle track with misrepair models. A limited number of ions are currently available, for protons, helium, carbon and nitrogen ions. Note Medras does not currently sanity-check input LET values, and will attempt to extrapolate based on available data even if unphysical tracks are specified (e.g. 500 keV/μm proton tracks). 

As an example, simulation of normal fibroblasts exposed to 80 keV/μm carbon ions can be simulated as:

```py
from medras import medrascell
testCell = medrascell.singleCell() # Initialise default cell object
exposure = {'dose':0, 'time':0, 'LET':80, 'particle':6} # Set up exposure
for d in [0,1,2,3,4,5,6]:
	exposure['dose'] = d # Can reuse exposure dictionary by updating dose value
	print(d,'\t',testCell.survival(exposure)) # Print survival for each dose
```

Finally, as noted above dose can be specified either as a single instantaneous exposure, or as a prolonged exposure. For prolonged exposures, each exposure (dose, time of exposure) is recorded in a nested list, with dose in Gy and time in hours. Thus, for example, to simulate a dose of 2 Gy delivered in 3 minutes, followed by a one-hour gap, followed by another dose of 2 Gy in 5 minutes, dose should be specified as: `[ [2,0.05],[0,1.0],[2,0.05] ]`. An illustration of the generation of a split-dose recovery curve is presented below:

```py
from medras import medrascell
testCell = medrascell.singleCell() # Initialise cell object
for t in [0,0.25,0.5,0.75,1,1.5,2,4,6]: # Time array specified in hours
	exposure = {'dose':[ [2,0.05],[0,t],[2,0.05]]} # 2 Gy in 5 minutes, plus variable gap
	print(t,'\t',testCell.survival(exposure)) 
```

### Other endpoints

Other non-survival endpoints can be modelled in a similar way, by calling the `modelExposure` method. This will return an list consisting of: 

- The current number of unrepaired breaks; 
- The current number of mis-repaired breaks;
- The current number of lethal chromosome aberrations;
- The current number of visible chromosome aberrations;
- The current number of mutations in the target gene; and
- The total number of DSB induced by the radiation exposure.

These values depend on both the dose delivered, and the time at which the endpoint is measured. As an example:

```py
from medras import medrascell
# Dose edependence of endpoints at 4 hours
testCell = medrascell.singleCell() # Initialise cell object
for d in [0,1,2,3,4,5,6]:
	exposure = {'dose':d,'time':2}
	print(d,'\t',testCell.modelExposure(exposure)) # Print endpoints as a function of dose

# Time dependence of endpoints at 2 Gy
for t in [0,0.25,0.5,0.75,1,1.5,2,4,6,12,24]:
	exposure = {'dose':2,'time':t}
	print(t,'\t',testCell.modelExposure(exposure)) # Print endpoints as a function of time
```

These endpoints are returned as a list. In addition, a `getFociCount` method can be used to calculate the number of foci at a given timepoint, taking into account the delay between DNA damage rejoining and foci clearance. 

### Model parameters

By default, the cell model is initialised with a set of model parameters which have been determined by fitting to experimental data (see (McMahon 2016, McMahon 2017) ). These parameters can also be varied when cell objects created, by passing in dnaParams and survParams dictionaries as arguments. 

The syntax of these objects can be seen in the defaults in `medrascell.py`, and their interpretation can be found in the referenced papers.

### Examples

Most of the above examples of different cell simulations are illustrated by the `medraasexamples.py` script in the `examples` folder, which simulates basic dose response curves for different cell cycle phases, irradiation with carbon ions at different LETs, and the kinetics of DNA repair, and chromosome aberrations following exposure to a range of doses, which should be possible to run directly from this folder. 

## Requirements

This code is written in python3, and requires the following libraries:

- numpy
- scipy

## Contacts

For questions/comments/bug reports, please contact stephen.mcmahon (at) qub.ac.uk

## References

[McMahon2017]	McMahon, S. J., McNamara, A. L., Schuemann, J., Paganetti, H., & Prise, K. M. (2017). A general mechanistic model enables predictions of the biological effectiveness of different qualities of radiation. Scientific Reports, 7(1), 688. http://doi.org/10.1038/s41598-017-10820-1

[McMahon2016]	McMahon, S. J., Schuemann, J., Paganetti, H., & Prise, K. M. (2016). Mechanistic Modelling of DNA Repair and Cellular Survival Following Radiation-Induced DNA Damage. Scientific Reports, 6, 33290. http://doi.org/10.1038/srep33290
