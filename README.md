# ConstraintSolverCMCBinaryKerr
Pseudospectral constraintsolver for binary head on Kerr-Black-Holes on CMC-Slices. This solves the setups discussed in https://arxiv.org/abs/0907.3163 and https://www.db-thueringen.de/receive/dbt_mods_00027284.

# Compilation

Execute the _Makefile_ with 

> make all

For advanced usage, in  _src/ftype.h_ the used datatype can be set e.g. _double_, _long double_. Afterwards the code has to be recompiled. 

# ConstraintSolver
The setup can be configured in _run/Config_.

- **Newton_itmin** minimal number of Newton-Iteration while solving the differential equation
- **Newton_itmax** maximum number of Newton-Iteration while solving the differential equation
- **Newton_tol_factor_steps** currently unused
- **Newton_tol_finish** acceptance tolerance for the Newton-Raphson-Method
- **Newton_verb** verbosity of the Newton-Raphson-Method, 1 for verbose, 0 for silent

- **num_threads** number of threads

- **ns_domain_0** number of grindpoints in radial direction in domain 0 (bordering the black hole)
- **ns_domain_1** number of grindpoints in radial direction in domain 1 (bordering the scri+)

- **nt** number of gridpoints in angular direction

- **Number_of_Solutions_in_the_Sequence** number of solutions in beetween the initial configuration and the final configuration

- **a0** distance of the black holes			
- **Rp** distance of scri+
- **Ktr** extrinsic curvature of the slice

- **r1** radius of the upper black hole, a marginal trapped surface is enforced here
- **c1** final free parameter of the upper black hole
- **Sz1** final spin of the upper black hole
- **Pz1** final momentum of the upper black hole

- **r2** radius of the lower black hole, a marginal trapped surface is enforced here
- **c2** final free parameter of the lower black hole
- **Sz2** final spin of the upper lower hole
- **Pz2** final momentum of the lower black hole

The code is run from within _run/_ with

> ../src/CMC_Data

The solver walks from c=Sz=Pz=0 to the configured values in _config/Config_. For every constructed solution the apparent horizon is searched. If found the physical properties of the horizon(s) are calculated.
The physical properties are stored in _run/physical_data. The data for constructing plots is found in _run/plots_. The raw data is found in _run/solution/. The different solutions of the constraint equation are identified by a timestamp which is attached to the filenames.
