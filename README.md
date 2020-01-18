# Companion repository for the paper "Fluid Approximation for Broadcasting Systems".

This repository contains script and classes needed to replicate the experiments described in the paper

Fluid Approximation for Broadcasting Systems
Luca Bortolussi, Jane Hillston and Michele Loreti

To replicate experiments the following tools are needed:
- JDK (to build and run simulation classes)
- Python 3 and SciPy (to plot graphs and solve ODE systems)

Once you have downloaded and installed the packages above you have to:


1. Clone git repository:

`git clone <url>`

2. Build java classes with gradle:

`./gradlew build`

3. Run simulations (this task may require many time):

`./gradlew run`

4. Run plotting scripts in folder `plotscript`:

`python plot_rb_files.py`

and

`python plot_gossip.py`




