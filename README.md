# ENGLISH
## Phoenix's Interior Ballistics Solver (PIBS)
A series of Python scripts (with a tkinter Graphical User Interface) for solving the interior ballistics system-of-equations, with provisions for both constrained design and certain optimization calculations.    

## Features
### Calculation from Parameters (Forward Calculation) 
* 0 dimensional interior ballistic SoE in the M.E. Serebryakov formulation.
    
    The interior ballistic problem is formulated in the orhtodox, traditional method as is taught in the Soviet Union and in China. The calculation is done in the reduced form for conventional and recoiless guns, and the expanded form for high-low launcher. 

* Efficient Integration

    Integrate up to user specified precision using high-order, adaptive [Runge Kutta Fehlberg](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method) 7(8)-th order method.

* Supports Shaped Propellant

    Baked shapes includes sphere, strip, cylinder, and multi-perforated cylinder/prism/rosette.

* Pre-defined Propellant Presets

    Baked in propellant definition for many propellants, sourced from the public literature, with power law, in the form of `a*p^b` burn rate-pressure relation, see `data/propellant.csv`. The burn rate is provided with a fudge factor to allow for granular modification to match experimental result, and can be used to simulate firing at off-normal temperature conditions, if the coefficients are known.


### Constrained Design/Optimization (Backward Calculation):
* Constrained Design conforming to Performance Characterstics

    The required web to achieve a desired peak pressure point can be solved either alone, with option `Lock Gun Length` in the GUI or passing `known_bore=True`, or with desired muzzle velocity, with option `Constrain Design`, in which case the required gun length will be solved simultaneousely, for both conventional and recoiless gun.

* Optimization conforming to Performance Characteristics

    From a desired system performance specification, the optimum loading fraction is searched and the one that minimizes total bore volume (the sum of chamber and bore volume) is found via varying the load fraction, with option `Minimize Tube Volume`, for both conventional and recoiless gun.

## How-To
* Use the calculator:
  - download the latest executable at [latest release](https://github.com/octo-org/octo-repo/releases/latest)
* Development setup:
  - Install Python (Developed and tested on 3.8.18, for Windows 7 compatibility)
  - Install packages via pip: matplotlib, matplotlib-label-lines
  - entry point is IB.py in /src

## Resources Used
* tcl/tk themes used include "awdark" & "awlight" from awthemes: https://wiki.tcl-lang.org/page/awthemes
* monospaced, CJK compatible font: Sarasa-Gothic, https://github.com/be5invis/Sarasa-Gothic
