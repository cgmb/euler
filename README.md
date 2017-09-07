```````````````````````````````````````````
 _|_|_|_|            _|
 _|        _|    _|  _|    _|_|    _|  _|_|
 _|_|_|    _|    _|  _|  _|_|_|_|  _|_|
 _|        _|    _|  _|  _|        _|
 _|_|_|_|    _|_|_|  _|    _|_|_|  _|

```````````````````````````````````````````

An ASCII art fluid simulator

```
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
X                                                              X
X                                                              X
X                                                              X
X                                                              X
X                                                              X
X                                                              X
X                                                              X
X                                                              X
X                                                              X
X               oO0O000                                        X
X           o0O0000000000                                      X
X       Oo0o00000000000O000                                    X
Xo o0O00000O00O000000000o000                                   X
X00000000o000000O000000000000O0o0OO                            X
X000000O00000000000o0000000000O0o0000OOo0oo                    X
X00000000000000000000000000000000o000O0000000O00oooo          oX
X00000000000000000000000000000O000000o0000000000000000000000000X
X00000000000000000000000000000O0000OO00000000000000000000000000X
X0000000000000000O0000000o000000000OOO0O000000000000O00000000O0X
XO0000000000000O0O0000000000000O0000O0o000000000000000000000000X
X000O0O00O0O0O00000000000000000000000000000O00000000000000O0000X
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
```

## About
Euler is a basic Eulerian fluid simulator based on the first edition of "Fluid
Simulation for Computer Graphics" by Robert Bridson. The book is an improved
version of his [SIGGRAPH 2007 Course Notes][1], which were developed in
conjunction with Matthias Müller-Fischer.

The project aim is to be minimal, correct and fun. It runs on Linux and OSX,
but unfortunately not on Windows. Sorry! If you're not able to run Euler, the
next best alternative is watching [a replay][ex1].

## How to Build

    cmake -H. -Bbuild && cmake --build build

## How to Run
Just pass a scenario file as an argument. Scenario files are plain text files,
where `0` represents fluid, `X` represents a solid wall, `?` represents a fluid
source and `=` represents a fluid sink. Maybe check out one of the
[sample scenarios](scenarios) for ideas.

    build/euler scenarios/block.txt

Don't tell anyone, but I heard there's also a secret `--rainbow` flag.
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
* Bridson recommends using RK2 integration for advection but this simulation
only uses RK1 integration, which necessitates smaller timesteps.

[1]: https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf
[2]: https://github.com/rlguy/GridFluidSim3D
[3]: https://github.com/austinEng/WebGL-PIC-FLIP-Fluid
[ex1]: https://asciinema.org/a/125371
[ex2]: https://asciinema.org/a/125380
