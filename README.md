# Ball Balancing Robot Control Project
This project explores control strategies for a ball balancing robot, focusing on two primary methods: Linear Quadratic Regulator (LQR) and Root Locus (PID).

## System Modeling and Equations of Motion
This section includes files related to the derivation of the system's equations of motion and their state-space formulation.
- `equations_of_motion.ipynb`: Jupyter notebook detailing the derivation of the system's equations of motion.

## Linear Quadratic Regulator (LQR)

The LQR method is used for optimal control design and state estimation in this project.
Files related to LQR:
- `ME4010_Project report.pdf`: Contains the theoretical background and numerical results for LQR controller design, observability analysis, and full-order/minimum-order Kalman observers.
- `LQR_ctrl.m`: MATLAB script for LQR controller implementation.
- `LQR_ctrl.slx`: Simulink model for the LQR controller.
- `LQR_obs.slx`: Simulink model for the LQR controller with full-order observer.
- `LQR_obs_noise.slx`: Simulink model for the LQR controller with full-order observer and measurement noise, filtering. 
- `LQR_min_obs.slx`: Simulink model for the LQR controller with minimum-order observer.
- `LQR_min_obs_noise.slx`: Simulink model for the LQR controller with minimum-order observer and noise filtering.
- `noise_data.mat`: Noise data captured for FFT.

## Root Locus Method (PID)

The Root Locus method is applied for designing classical cascade PID control architecture.
Files related to Root Locus:
- `ME4010_Project report.pdf`: Details the cascade control strategy and PID controller tuning using Root Locus analysis.
- `Rlocus_ctrl.m`: MATLAB script for Root Locus based PID controller implementation.
- `Rlocus.slx`: Simulink model for the Root Locus based PID controller.
