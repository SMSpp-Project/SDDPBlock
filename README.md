# SDDPBlock

SMS++ `:Block` and `Solver` for multistage stochastic programming problems
amenable of solution by means of Stochastic Dual Dynamic Programming (SDDP)
approaches. In particular three components are provided:

- `SDDPBlock` derives from `Block` and represents a multistage stochastic
  programming problem amenable of solution by means of SDDP approaches. This
  roughly speaking means that each stage only depends on decisions from the
  previous one. An `SDDPBlock` has a time horizon *T*, and *T* sub-`Blocks`,
  each one being a [StochasticBlock](https://gitlab.com/smspp/stochasticblock).
  The t-th sub-`Block` represents an approximation to the problem associated
  with stage t, in the sense that each stage (save for the last) must have
  already inside it a `PolyhedralFunction` that represents the currently
  available approximation of the future (average) cost of all the subsequent
  stages. `SDDPBlock` must be provided to pointers to these *T* - 1
  `PolyhedralFunction`. The stochastic nature of the problem is represented
  by a (currently, finite) set of *scenarios*, each if which is a vector of
  double and spans all the time horizon *T*. Each vector is divided into *T*
  parts, each one being associated with a stage of the multistage problem.

- `SDDPSolver` derives from `Solver` and implements the SDDP method for
  multistage linear stochastic problems by acting as a wrapper for the SDDP
  solver implemented in the open-source [STochastic OPTimization library
  (StOpt)](https://gitlab.com/stochastic-control/StOpt). `SDDPSolver`
  requires the problem to be specified as a `SDDPBlock`. The SDDP approach
  is basically a nested Benders' one, whereby the stages are solved in
  order 1, 2, ..., *T* with the currently available approximation of the
  future (average) cost of all the subsequent stages provided by the
  current state of the `PolyhedralFunction` in the `SDDPBlock`, for a
  properly chosen subset of the scenarios. The dual variables of the
  constraints in stage *t* + 1 containing the decisions of stage *t* are
  then used to compute a valid linear lower approximation of the cost of
  stage *t* + 1 (that includes the cost of all subsequent stages, if any)
  which are (averaged over the scenarios and) added to the function of
  stage *t*. These are basicallt Benders' cuts, in fact the process is
  realised via a `BendersBFunction`. The approach iterates these *forward
  passes* (taking decisions based on the current functions) and *backward
  passes* (updating the functions bases on the outcome of previous decisions)
  until convergenge is attained. As in Benders' decomposition this requires
  dual variables, and hence the subproblems corresponding to the stages to be
  convex (or some sort of convex relaxations of them being solved).

- `SDDPGreedySolver` derives from `Solver` and implements a sequential,
  greedy strategy to solve an `SDDPBlock` for a fixed scenario. Its typical
  use is to evaluate the decisions taken by `SDDPSolver` by providing a
  bunch of feasible (although not necessarily optimal) solutions for each
  scenario; doing this for many scenarios (possibly some not considered
  by `SDDPSolver`) provides an upper bound on the expected cost of the
  problem. The idea of `SDDPGreedySolver` is to avoid to solve the
  whole multi-stage problem (that can be very large even for one fixed
  scenario if *T* is large and/or the `:Block` corresponding to each stage
  are) but rather exploit the future (average) cost of all the subsequent
  stages stored in the `PolyhedralFunction` of the `SDDPBlock`. Basically
  the problem corresponding to the (`:Block` representing) each stage is
  solved in isolation, for fixed values of the decisions in the previous
  stage (if any), assuming that the future (average) cost function of all
  the subsequent stages provides appropriate guidance for the solution of
  the current stage without a need of explicitly looking at the detailed
  future decisions. This is typically done *after* that `SDDPSolver` has
  been ran on the `SDDPBlock`, so that the crucial future (average) cost
  functions are hopefully "good enough". Yet, the process is clearly an
  heuristic one and does not exactly the deterministic (single-scenario)
  multistage problem, even less the stochastic problem encoded by the
  `SDDPBlock`.


## Getting started

These instructions will let you build `SDDPBlock` on your system.

### Requirements

- [StochasticBlock](https://gitlab.com/smspp/stochasticblock),
  which in turn requires the "core SMS++"

- [STochastic OPTimization library (StOpt)](https://gitlab.com/stochastic-control/StOpt)

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(SDDPBlock)
target_link_libraries(<my_target> SMS++::SDDPBlock)
```

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `SDDPBlock` has to include a "main makefile" of the
module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). Relevant examples are the
[sddp_solver](https://gitlab.com/smspp/tools/-/blob/develop/sddp_solver/sddp_solver.cpp?ref_type=heads)
available in the [tools](https://gitlab.com/smspp/tools) repository.
The makefiles in turn recursively include all the required other makefiles,
hence one should only need to edit the "main makefile" for compilation type
(C++ compiler and its options) and it all should be good to go. In case some
of the external libraries (say, StOpt) are not at their default location, it
should only be necessary to create the `../extlib/makefile-paths` out of the
`extlib/makefile-default-paths-*` for your OS `*` and edit the relevant bits
(commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.


## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/sddpblock/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Rafael Durbano Lobato**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
