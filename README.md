# PhoenixYard
A grpahical user interface for solving the internal ballistics system-of-equation made using tkinter.

# Feature
* Forward Calculation:
  - Internal ballistic SoE in the M.E. Serebryakov formulation.
  - Integrate up to user specified precision using adaptive Runge Kutta Fehlberg 7(8)-th order method.
  - Propellant shape: sphere, strip, cylinder, multi-perforated cylinder/prism/rosette.
  - Supports custom propellant definition with a x p^b burn rate-pressure relation. (Currently the propellant csv file is wrapped into the binary in the one-file mode, therefore its necessary to rebuild the program to incorporate updated data. Work is on going.) 
* Constrained Design/Optimization:
  - Solve the web thickness and length of gun requried to achieve a peak pressure & muzzle velocity specification
  - Solve the web thickness, length of gun and load factor required to achieve a peak pressure & muzzle velocity specification while minimizing chamber volume.

# How-To
* Compiled executables:
  - clone the repository and use the compiled windows executables in /bin
* Development setup:
  - Install Python (>3.8)
  - Install packages via pip: matplotlib, matplotlib-label-lines
  - entry point is app.py in /src

# Resources Used
* tcl theme used "awdark" & "awlight" from awthemes: https://wiki.tcl-lang.org/page/awthemes
* Hack font, from https://github.com/source-foundry/Hack
