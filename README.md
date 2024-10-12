# Simple Pack Modeling Interface w/PyBAMM

Lithium-ion batteries have been ubiquitious for numerous energy applications. In particular, there is great interest in researching battery packs (i.e. simulating modules with multiple cells in series/parallel). Often, reduced-order models are preferred for pack-level studies so as to reduce computational overhead.

Along these lines, this repository seeks to accomplish the following: 
1) Detail an easy-to-use, robust API for simulating 'variable-size' battery packs 
2) Allow modelers to rapidly conduct cycling studies with variated operating conditions and parameters
3) Provide a post-processing and plotting interface to readily wrangle with large volumes of simulation results
4) Serve as a tool for research groups at large
    - _Developed primarily for INSERT OUR GROUP AND PAPER LINK HERE_ 
5) Serve as a well-documented tutorial/exercise for simulation development
    - _Existing libraries, though supportive of complex models, can be difficult to comprehend for beginner programmers and electrochemists alike. This codebase strives to be a simple introduction_

**IMPORTANT NOTE:** THIS CODEBASE IS VERY MUCH A WORK-IN-PROGRESS. ALTHOUGH FUNCTIONAL, THERE IS CERTAINLY ROOM FOR STRUCTURAL IMPROVEMENTS 


## Configuring Simulations
`mainmodel.py` provides an example of a pack simulation setup

### Operating Conditions
The following **operating conditions** are enumerated at the top of `mainmodel.py`

| **Parameter**           | **Description**                                                                 | **Example**                                |
|-------------------------|---------------------------------------------------------------------------------|------------------------------------------|
| NUM_SERIES              | Number of cells in SERIES                                                       | 5                                      |
| NUM_PARALLEL            | Number of cells in PARALLEL                                                     | 5                                      |
| NUM_CYCLES              | Number of cycles (full discharge + CC-charge + CV-charge is one cycle)          | 100                                    |
| C_RATE                  | Rate of CC-discharge/charge. Applied current computed from this                 | 1.0                                   |
| VOLTAGE_WINDOW         | Low and High Cutoff voltages                                      | (2.5 * NUM_SERIES, 4.1 * NUM_SERIES)                         |
| CURRENT_CUT_FACTOR      | Stop condition for CV-charge; if the solved pack current <= CUT_FACTOR*I_INPUT, halt! | 1/10                                     |
| CAPACITY_CUT_FACTOR      | Stop condition for capacity fade. If ANY cell's capacity is CUT_FACTOR of original capacity, halt! | 0.8                                     |
| DISCRETE_PTS            | How many points in particle mesh                                                | 30                                       |
| HOURS                   | Duration of simulation. Ideally, derived from C-rate                            | 2                                        |
| TIME_PTS                | Number of time points to return solution PER charge/discharge                   | 100                                      |
| EXPERIMENT              | Name of study. Each study should get a unique name; all data outputted to namesake folder | "5by5_100cycles_const"          |

`USE_C_RATE = True`,  `C_RATE` value is used to compute **applied pack current**
`USE_C_RATE = False`, `I_INPUT` value is used AS the **applied pack current**


### Parameter Variation
Model parameters can be varied across all the cells in a pack. This is useful for characterizing performance of heterogenous battery modules.

The following parameters (for simple SPM model) are enumerated in `params.py`

| **Parameter**           | **Description**              | **Example Value** |
|-------------------------|------------------------------|-------------------|
| **POS_DIFFUSION**       | Cathode Diffusion Coefficient            | 1.0e-14           |
| **POS_CSN_MAX**         | Cathode Max Concentration     | 51555             |
| **POS_CSN_INITIAL**     | Cathode Initial SOC (mol/m3)        | 51555*0.5         |
| **POS_ELEC_THICKNESS**  | Cathode Thickness (m)          | 80e-6             |
| **POS_ELEC_POROSITY**   | Cathode Porosity  (%)           | 0.385             |
| **NEG_DIFFUSION**       | Anode Diffusion Coefficient              | 2.0e-14           |
| **NEG_CSN_MAX**         | Anode Max Concentration      | 30555             |
| **NEG_CSN_INITIAL**     | Anode Initial SOC (mol/m3)           | 30555*0.74        |
| **NEG_ELEC_THICKNESS**  | Anode Thickness (m)             | 88e-6             |
| **NEG_ELEC_POROSITY**   | Anode Porosity  (%)             | 0.485             |
| **PARTICLE_RADIUS**     | Particle Radius (m)             | 2e-06             |

_For further remarks on how parameters 'fit' into the model equations/DAE, see the `Electrochemical Model POV` section below_
_For source of parameter values, see the `References` section below_ **include pointer to paper**
 
Variations can be applied using a few different schemes. See the following examples with cathode porosity:

**Notations**  
$\epsilon_k$ = porosity of a given cell 'k'  
$v$ = porosity value  
$p$ = porosity percent  
$\sigma$ = porosity standard deviation

| **Code**                                                     | **Distribution Notation**                                        |
|--------------------------------------------------------------|--------------------------------------------------|
| `POROSITY = Variator.from_percent("", v, p)` | $\epsilon_k \sim Uniform(v - \frac{p\dot v}{100},v+\frac{p\dot v}{100})$ |
| `POROSITY = Variator.from_gaussian_percent("", v, p)` | $\epsilon_k \sim N(v, \frac{p\dot v}{100})$          |
| `POROSITY = Variator.from_gaussian_stddev("", v, σ)` | $\epsilon_k \sim N(v,σ)$                   |

_The first argument (empty string) is an arbitrary name that can be given to the parameter (no functional importance)_

## Data Output
Each `experiment` is outputted to namesake folder under `data/`.  
`data/EXAMPLE/` provides an example of a simulation study output (all files generated from a SINGLE experiment)

| Name          | Info                                                                                                 |
|---------------|------------------------------------------------------------------------------------------------------|
| capacities.csv| Discharge capacity of EACH cell after EACH discharge cycle                                         |
| data.csv      | Master simulation data. <br> -Cols 1-3: Cycle #, Protocol, Time Index. <br> -Cols for **top-level** attributes: Pack Voltage, Pack Current, String Currents ('String' is a chain of cells in series). <br> -Cols for **cell-level** attributes: Concentration SOC, SEI Length, Voltage, Capacity Integration. |
| profile.json  | Simulation attributes, operating conditions, applied parameter variations enumerated                 |
| model.pkl     | The "Pack" object (src/pack.py). Pickled/unpickled to access internal attributes                      |

**Cell naming convention:**  
A cell in **3rd** 'string' in parallel and **2nd** cell in series: **Cell 2,3**

## Post-processing and Plotting
`reader.py` provides an example of the data post-processing interface  
_Refer to comments in file until further documentation written... TBD_

## Developers' Guide

### Electrochemical Model POV 

#### Single-Particle-Model with SEI

### Programmer POV

## References


