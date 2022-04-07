# SC42125
## By Frédéric Larocque

## What is this?
Matlab code that supports the modeling and design of MPC controllers for the longitudinal control of a H145 helicopter. The 3DoF longitudinal model of the helicopter was derived and linearized around hover. Three types of MPC controllers were designed:
- Regulation controller
- Trajectory tracking controller (state feedback)
- Offset-free output tracking controller

## How to use?
Run each ".m" file to run a specific MPC controller.

## File list
- functions : folder containing the functions supporting all the main files
- figures: folder containing the figures used in the report
- papers : folder containing the documents, lecture notes and presentations used to complete the assignment
- compare_linear_non_linear.m: compares the dynamics of the non-linear and linearized model around the hover point
- init_hover_dynamics.m: file to initialize the dynamics of the helicopter
- init_hover_dynamics_robustness_test.m: same file as init_hover_dynamics.m, but with additionnal mass
- offset_free_tracking.m: offset-free tracking MPC controller
- reference_controller.m: trajectory tracking MPC controller
- regulation_MPC.m: regulation MPC controller and its comparison to the LQR controller
- robustness.m: not used in report, but analyzes behaviour of MPC when dynamics are not exactly as the MPC's internal model
- stability_analysis.m: analyzes stability, attraction set and terminal set of MPC controller
