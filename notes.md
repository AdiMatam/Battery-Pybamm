# Setting up a Model -- High Level Concept:

1. Initialize a model (BaseModel is fundamental starting point) 
2. Define variables to observe and parameters 
    - pybamm.Variable() --> to be computed
3. Define equations for each variable -- passed into `model.rhs` dictionary
4. Define `initial_conditions`, `boundary_conditions`, `variables` for `model`
5. Process model, geometry (when applicable) with parameter inputs
6. Create mesh and discretize geometry
7. Run solver
    - Solution evaluates each variable @ each time step (specificed during mesh discretization) @ each point in spatial dimensions
    - solution.y --> Array<Array<int, TimeSteps>, Pts>


# Thoughts

- May want wrappers for some of these... 
    -> Boundary Conditions
    -> Geometry
    -> Too many string based parameters / variables. Difficult to maintain

# Done @ 08/14
- Rough meeting notes on SPM parameters (porosity etc) and capacity calculation 
- Fix/Debug model --> removed polar term (1/r2) as already accounted for in "spherical" domain paramter

# Done @ 08/06

- Study PDE script
    - How does parameterization work? Radius? [Here](https://docs.pybamm.org/en/latest/source/examples/notebooks/parameterization/parameterization.html)

- Start creating SPM model
- Types of `submeshes` and `spatial methods`

- Building models w/out guidance -- going from mathematical notation -> pybamm solution
    -> pde and ode examples. with pde, figured out how paramaterization processing is done

- Began building SPM model
    -> go over ## items

- Interpreting the real equation(s) -- boundary conditions in particular! 