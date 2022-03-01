# Fretting dynamics
Simulating various friction systems. Note that this is very much a work in progress; use on own risk.

Current features: Simulates 1-dimensional friction chains, like Burridge-Knopoff systems. Calculates friction evolution and "noise" to pressure and see evolution of partial-slip. Various simplistic tools from chaos theory to study state space orbits etc.

## TODO
- [ ] Write test cases (fuzzing?)
- [X] 1DOF rate-dependent friction
  - [X] Rate-dependent friction
  - [X] Poincaré maps
  - [X] Frequency analysis
- [X] N-block system
- [ ] Refactor (see design philosophy further down)
- [ ] Write benchmarks
- [ ] YAML-parser
- [ ] Use fmt-library instead of native string concat
- [ ] Simulation types
  - [X] Simple time integrator
  - [ ] Frequency sweep
  - [ ] Snapshots throughout a long sim
- [ ] HertzRateSine
  - [X] Poincare sections
  - [X] Slip history and FFT
  - [ ] Roughness (in pressure), e.g. using Ra/Rz
- [ ] Monte-Carlo simulation
- [ ] Change to Conan-based build?

## Design Philosophies
- Const and constexpr/eval of as much as possible
- As little dynamic dispatching as possible
- Functionality before composability
- Use builder-pattern?
- For heap objects: single owner and unique_ptr


## Licence
Free or whatever? GPL
