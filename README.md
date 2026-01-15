# Description
This program walkthrough documents the trajectory optimization and mission design for an interplanetary intercept of the asteroid 99942 Apophis. This analysis identifies the most fuel-efficient launch and arrival windows using high-fidelity orbital data and a numerical solver.

## 1. Project Overview & Motivation

The objective was to apply Lambert’s problem-solving techniques to a trajectory from Earth to Apophis within a specific launch window from early April to early June 2028. This project moves beyond simple 2D orbit by building a MATLAB script capable of processing vast ephemeris data to find a global minimum for fuel consumption ($\Delta V$) and maximize scientific payload delivery.

## 2. Implementation of the Lambert Solver

The primary component of this mission design is the `Lambert_aro3090` function, which determines the unique orbit connecting two points in space for a specified time of flight (TOF).

### Mathematical Foundation
The function is a hybrid implementation of the Gooding (1990) and Lancaster-Blanchard (1969) procedures. It solves for the universal variable $x$, which represents the energy level of the orbit (elliptical, parabolic, or hyperbolic).

### Segmented Vectorization Logic
To handle the scale of the search space—iterating through thousands of potential launch and arrival dates at 0.1-day timesteps—the function utilizes segmented vectorization.
<br />

* **Memory Management:** The function splits input data into segments (controlled by the `nm` variable) to process millions of potential trajectories without exceeding memory limits.
* **Multi-Revolution Capability:** By using an imaginary number input for `Nrev`, the function automatically identifies all possible solutions (prograde and retrograde) for multiple revolutions around the Sun.

## 3. Methodology: Iterative Optimization

The trajectory design relies on a nested loop architecture to construct a comprehensive map of potential mission profiles.
<br />

* **Data Acquisition:** 3D position and velocity vectors for Earth and Apophis are loaded for the 2028–2029 period.
* **Velocity Calculation:** For every date combination, the Lambert function calculates the required transfer velocities ($V_1$ and $V_2$).
* **Cost Function Evaluation:** The total $\Delta V$ is calculated as the sum of the departure maneuver (the difference between the transfer $V_1$ and Earth's velocity) and the arrival maneuver (the difference between the transfer $V_2$ and the asteroid's velocity).

## 4. Results & Mission Characteristics

The simulation successfully converged on an optimal mission profile that maximizes dry mass delivery while meeting all timing constraints.

<p align="center">
  <img src="https://i.imgur.com/glMZS23.png" height="80%" width="80%" alt="Apophis Mission Trajectory Plot"/>
  <br />
  <i>Figure 1: Optimal Interplanetary Transfer Trajectory to Asteroid Apophis</i>
</p>

### Final Mission Deliverables
<br />

* **Optimal Launch/Arrival:** The solver identified the specific dates in 2028 and 2029 that provide the lowest energy path.
* **Characteristic Energy ($C_3$):** Calculated at departure to determine the wet mass capacity of the New Glenn launch vehicle based on its specific performance curves ($K_1, K_2, K_3$).
* **Mass Efficiency:** Using the Tsiolkovsky Rocket Equation, the code determines the final scientific dry mass and the required propellant mass based on the engine's specific impulse ($I_{sp}$).

By combining numerical integration with launch vehicle performance modeling, this project provides a complete end-to-end mission architecture for the intercept of Apophis.
