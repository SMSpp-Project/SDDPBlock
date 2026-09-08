# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added 

### Changed 

- `CutProcessing` is gone: removing the redundant rows of a
  `PolyhedralFunction` is done by the `PolyhedralFunctionBlock` that owns it,
  through `PolyhedralFunctionBlock::remove_redundant_rows()`

### Fixed 

## [0.5.4] - 2025-12-12

### Added 

- SDDPSolver parameters to better control outputs

- parameter to set the seed of the forward simulator

### Changed 

- adapted to new standard organization of makefiles

- parameter in set\_ComputeConfig() is now const

- made USE_MPI public

### Fixed 

- corrected dysfunctional makefile-s and -c

## [0.5.3] - 2024-02-29

### Changed 

- adapted to new CMake / makefile organisation

## [0.5.2] - 2023-17-05

### Added

- Handle eEverykIteration events in SDDPSolver.

### Fixed

- Windows compilation issues.

## [0.5.1] - 2022-07-01

### Added

- SDDPGreedySolver::set\_ComputeConfig()

- Load cuts from file

- Allow scenarios to be randomly chosen in SDDPGreedySolver

- Multiple parameters to SDDPGreedySolver

### Changed

- Update interface with StOpt

- SDDPGreedySolver becomes a CDASolver

- Output of SDDPSolver and SDDPGreedySolver

- Define the sense of the "Objective" of the SDDPBlock

## [0.5.0] - 2021-12-08

### Added

- Multiple parameters to SDDPSolver and SDDPGreedySolver

- Support for multiple sub-Blocks per stage in SDDPBlock

- ParallelSDDPSolver

- SDDPSolverState

- SDDPSolver::set_ComputeConfig()

- Handling Configuration for get\_var\_solution()

- Storage of random cuts

### Changed

- Cuts provided by StOpt are added incrementally

### Fixed

- Objective value of subproblem in SDDPGreedySolver

- Bug in oneStepForward() regarding the simulation id

## [0.4.0] - 2021-05-02

### Added

- Inner Blocks of SDDPBlock can be configured by SDDPSolver and
  SDDPGreedySolver

- Implementation of SDDPSolver

## [0.3.0] - 2020-09-16

### Added

- SDDPGreedySolver

## [0.2.0] - 2020-03-06

### Added

- Serialization

## [0.1.0] - 2020-01-09

### Added

- First test release

[Unreleased]: https://gitlab.com/smspp/sddpblock/-/compare/0.5.4...develop
[0.5.4]: https://gitlab.com/smspp/sddpblock/-/compare/0.5.3...0.5.4
[0.5.3]: https://gitlab.com/smspp/sddpblock/-/compare/0.5.2...0.5.3
[0.5.2]: https://gitlab.com/smspp/sddpblock/-/compare/0.5.1...0.5.2
[0.5.1]: https://gitlab.com/smspp/sddpblock/-/compare/0.5.0...0.5.1
[0.5.0]: https://gitlab.com/smspp/sddpblock/-/compare/0.4.0...0.5.0
[0.4.0]: https://gitlab.com/smspp/sddpblock/-/compare/0.3.0...0.4.0
[0.3.0]: https://gitlab.com/smspp/sddpblock/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/sddpblock/-/compare/0.1.0...0.2.0
[0.1.0]: https://gitlab.com/smspp/sddpblock/-/tags/0.1.0
