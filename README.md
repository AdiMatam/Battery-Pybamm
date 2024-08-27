# Simple Pack Modeling Interface

Lithium-ion batteries have been ubiquitious for numerous energy applications. In particular, there is great interest in researching battery packs (i.e. simulating modules with multiple cells in series/parallel). Often, reduced-order models are preferred for pack-level investigation, so as to reduce computational overhead.

Along the aforementioned lines, the following repository seeks to accomplish several targets:
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
| THEORETICAL_CAPACITY     | Theoretical capacity per unit area in A/m² of the cells                         |27.263836618115                    |
| NUM_SERIES              | Number of cells in SERIES                                                       | 5                                      |
| NUM_PARALLEL            | Number of cells in PARALLEL                                                     | 5                                      |
| NUM_CYCLES              | Number of cycles (full discharge + CC-charge + CV-charge is one cycle)          | 100                                    |
| C_RATE                  | Rate of CC-discharge/charge. Applied current computed from this                 | 1.0                                   |
| BASE_CURRENT            | Estimated capacity * C_RATE                                                     | _Derived_                                      |
| I_INPUT                 | BASE_CURRENT * NUM_PARALLEL                                                     | _Derived_                                      |
| CURRENT_CUT_FACTOR      | Stop condition for CV-charge; if the solved pack current <= CUT_FACTOR*I_INPUT, halt! | 1/10                                     |
| VOLTAGE_LOW_CUT         | Cutoff voltage at which to stop discharging                                      | 2.5 * NUM_SERIES                         |
| VOLTAGE_HIGH_CUT        | Cutoff voltage at which to stop CC-charging                                      | 4.1 * NUM_SERIES                         |
| DISCRETE_PTS            | How many points in particle mesh                                                | 30                                       |
| HOURS                   | Duration of simulation. Ideally, derived from C-rate                            | 2                                        |
| TIME_PTS                | Number of time points to return solution PER charge/discharge                   | 100                                      |
| EXPERIMENT              | Name of study. Each study should get a unique name; all data outputted to namesake folder | "5by5_100cycles_const"          |


### Parameter Variation
Model parameters can be varied across all the cells in a pack. This is useful for characterizing performance of heterogenous battery modules.

The following parameters (for simple SPM model) are enumerated in `params.py`

| **Parameter**           | **Description**              | **Example Value** |
|-------------------------|------------------------------|-------------------|
| **POS_DIFFUSION**       | Cathode Diffusion            | 1.0e-14           |
| **POS_CSN_MAX**         | Cathode Max Concentration     | 51555             |
| **POS_CSN_INITIAL**     | Cathode Initial SOC          | 51555*0.5         |
| **POS_ELEC_THICKNESS**  | Cathode Thickness            | 80e-6             |
| **POS_ELEC_POROSITY**   | Cathode Porosity             | 0.385             |
| **NEG_DIFFUSION**       | Anode Diffusion              | 2.0e-14           |
| **NEG_CSN_MAX**         | Anode Max Concentration      | 30555             |
| **NEG_CSN_INITIAL**     | Anode Initial SOC            | 30555*0.74        |
| **NEG_ELEC_THICKNESS**  | Anode Thickness              | 88e-6             |
| **NEG_ELEC_POROSITY**   | Anode Porosity               | 0.485             |
| **PARTICLE_RADIUS**     | Particle Radius              | 2e-06             |

_For further remarks on how parameters 'fit' into the model equations/DAE, see the `Electrochemical Model POV` section below_
 
Variations can be applied using a few different schemes. See the following examples with cathode porosity:

**Notations**  
$\epsilon_k$ = porosity of a given cell 'k'  
$v$ = porosity value  
$p$ = porosity percent  
$\sigma$ = porosity standard deviation

| **Code**                                                     | **Distribution Notation**                                        |
|--------------------------------------------------------------|--------------------------------------------------|
| `POS_ELEC_POROSITY = Variator.from_gaussian_percent("any name", v, p)` | $\epsilon_k \sim Uniform(v - \frac{p*v}{100}, v + \frac{p*v}{100})$ |
| `POS_ELEC_POROSITY = Variator.from_gaussian_percent("any name", v, p)` | $\epsilon_k \sim N(v, \frac{p*v}{100})$          |
| `POS_ELEC_POROSITY = Variator.from_gaussian_stddev("any name", v, \sigma)` | $\epsilon_k \sim N(v, \sigma)$                   |


## Post-processing and Plotting


## Developers' Guide

### Electrochemical Model POV 

#### Single-Particle-Model with SEI

### Programmer POV

## References


