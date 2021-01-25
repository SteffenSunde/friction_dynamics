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
- [ ] 1DOF rate-dependent friction
  - [X] Rate-dependent friction
  - [ ] Poincaré maps
  - [X] Frequency analysis
- [ ] N-block system
- [ ] Refactor (see design philosophy further down)
- [ ] Write benchmarks
- [ ] YAML-parser
- [ ] Simulation types
  - [ ] Simple time integrator
  - [ ] Frequency sweep
  - [ ] Snapshots throughout a long sim
- [ ] Models
  - [ ] HertzRateSine
  - [ ] 

## Design Philosophies
- Constant (compile-time evaluation) of as much as possible
- As little dynamic dispatching as possible
- Functionality before composability
- Use builder-pattern?


## Licence
Free or whatever?