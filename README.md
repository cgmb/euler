```````````````````````````````````````````
 _|_|_|_|            _|
 _|        _|    _|  _|    _|_|    _|  _|_|
 _|_|_|    _|    _|  _|  _|_|_|_|  _|_|
 _|        _|    _|  _|  _|        _|
 _|_|_|_|    _|_|_|  _|    _|_|_|  _|

```````````````````````````````````````````

An ASCII art fluid simulator

[![block of water](https://asciinema.org/a/125371.png)][ex1]

## About
Euler is a basic Eulerian fluid simulator based on the first edition of "Fluid
Simulation for Computer Graphics" by Robert Bridson. The book is an improved
version of his [SIGGRAPH 2007 Course Notes][1], which were developed in
conjunction with Matthias Müller-Fischer.

The project aim is to be minimal, correct and fun. It runs on Linux and OSX,
but unfortunately not on Windows.

## How to Build

    cmake -H. -Bbuild && cmake --build build

## How to Run
Just pass a scenario file as an argument. Scenario files are plain text files,
where `0` represents fluid, `X` represents a solid wall, `?` represents a fluid
source and `=` represents a fluid sink. Maybe check out one of the
[sample scenarios](scenarios) for ideas.

    build/euler scenarios/block.txt

Don't tell anyone, but I heard there's also a `--rainbow` flag you can pass.
[I wonder what that does...][ex2]

## Other Simulators
There are some pretty cool simulators out there. [GridFluidSim3D][2] by Ryan
Guy and [WebGL-PIC-FLIP-Fluid][3] by Austin Eng appear to be based off the
second edition of Bridson's book. They are more featureful and much faster than
Euler due to GPU acceleration, but the underlying fluid simulation is similar.

## Bugs and Limitations
* Fluid volume is not conserved. This can be corrected for, but conservation
of volume, mass and momentum are not inherent properties of the simulation.
* This implementation uses marker particles to track the fluid. They are rather
noisy in comparison to level set methods and they tend to cluster.
* Bridson recommends using RK2 integration for advection, but this simulation
only uses RK1 integration. That necessitates smaller timesteps.
* No preconditioning is done for the conjugate-gradient pressure solve, so
the iteration limit is hit quite regularly.

[1]: https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf
[2]: https://github.com/rlguy/GridFluidSim3D
[3]: https://github.com/austinEng/WebGL-PIC-FLIP-Fluid
[ex1]: https://asciinema.org/a/125371
[ex2]: https://asciinema.org/a/125380
