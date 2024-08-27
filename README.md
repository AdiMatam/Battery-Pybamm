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

## Post-processing and Plotting


## Developers' Guide

### Electrochemical Model POV 

#### Single-Particle-Model with SEI

### Programmer POV

## References


