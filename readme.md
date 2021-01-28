# Fretting dynamics
Simulating various friction systems

Plan:
- FFT of total shear force
- Implement friction law
- Study evolution, find some case of shakedown
- Implement random noise to contact pressure based on wavenumber and Ra?
- Study contact evolution
- Perform Monte Carlo-simulation of friction to study shakedown

## TODO
- [ ] Write test cases
- [X] 1DOF rate-dependent friction
  - [X] Rate-dependent friction
  - [X] Poincaré maps
  - [X] Frequency analysis
- [X] N-block system
- [ ] Refactor (see design philosophy further down)
- [ ] Write benchmarks
- [ ] YAML-parser
- [ ] Simulation types
  - [X] Simple time integrator
  - [ ] Frequency sweep
  - [ ] Snapshots throughout a long sim
- [ ] HertzRateSine
  - [X] Poincare sections
  - [X] Slip history and FFT
  - [ ] Roughness (in pressure)
- [ ] Monte-Carlo simulation

## Design Philosophies
- Const and constexpr/eval of as much as possible
- As little dynamic dispatching as possible
- Functionality before composability
- Use builder-pattern?
- For heap objects: single owner and unique_ptr


## Licence
Free or whatever?