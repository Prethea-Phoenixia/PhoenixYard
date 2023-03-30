# PhoenixYard
A grpahical user interface for solving the internal ballistics system-of-equation made using tkinter.

# Feature (Current and Planned)
* Numerically (Runge Kutta Fehlberg 4-5 th order) integrates the internal ballistics exponential burn rate SoE.
  - currently the burn rate is modeled as b*p^n therefore 0 shot start pressure is not yet supported.
* (WIP) Traditional analytical solution for multi-perf propellant.
  - optimisation with constraints planned for future.
* Editable propellant definition.
  - currently the geometry of the propellants are not editable, although it is believed that the ones present covers most use cases.
  - support for simpler geometries (grain, non-perforated tubes) planned for the future.

# How to Use
custom propellants can be loaded from propellants.csv. Otherwise just run the program.

# Image
![image](https://user-images.githubusercontent.com/42470911/228233534-efe123ab-ca39-4da0-bc8c-eeacad118425.png)

The credit for the modern looking UI goes to the awthemes team. The font is Hack.
