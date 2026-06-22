# reactiondiffusion
MATLAB scripts for 1D stochastic reaction-diffusion simulations with active Ornstein-Uhlenbeck noise.

## Scripts

generate_data_1d(points, runtime, t_corr, runs)
- points: number of spatial grid points
- runtime: number of time steps
- t_corr: active-noise correlation time
- runs: number of independent simulations
- use runs == 1 to plot one simulation and detect patches
- use runs > 1 to run repeated simulations
- example: generate_data_1d(1000, 5000, 5, 1)

generate_data_1d_data_collect(points, runtime, t_corr, patternpoints)
- points: number of spatial grid points
- runtime: number of time steps per simulation
- t_corr: active-noise correlation time
- patternpoints: number of detected patches to collect
- repeats simulations until enough patches are found
- example: generate_data_1d_data_collect(1000, 5000, 5, 100)

periodic_test(points, runtime, t_corr, runs)
- like generate_data_1d, but yields correct size for patches wrapping around boundary, fixes issue with regionprops height
- filters out patches that get cutoff by runtime
- example: periodic_test(1000, 5000, 5, 1)

periodic_test_selftest()
- checks the circular-height and boundary-filter helpers

x_first_passage/x_first_passage_tcorr_sweep.py
- reduced X first-passage sweep with reaction and diffusion removed
- writes a table, density CSV, and overlay plot
- example: python3 x_first_passage/x_first_passage_tcorr_sweep.py

generate_data_1d is for inspecting individual runs or repeated mean trajectories.
Repeated runs in generate_data_1d collect the mean activator trajectory over time across independent stochastic runs, not patch statistics.
generate_data_1d_data_collect is for repeatedly sampling patch statistics.
periodic_test is for checking patch geometry under periodic spatial boundaries.

generate_data_1d writes xmat.csv and ymat.csv for single runs, or simdat1d_*.csv for repeated runs.
generate_data_1d_data_collect writes dataforhist.csv.
periodic_test writes periodic_test_xmat.csv, periodic_test_ymat.csv, periodic_test_data.csv, and periodic_test_summary.csv for single runs.

## Parameters
Edit model parameters inside rd_step_active in each file:

```matlab
DX = 1;              % activator diffusion
DY = 5;              % inhibitor diffusion
t_v = 20;            % tau_v
gamma = 1/t_v;
betavar = 0.7*gamma; % beta = 0.7
alphavar = 0.5*gamma;% alpha = 0.5
epsilon = 1;
a = 0.1*sqrt(epsilon);
```
Some scripts currently use different t_v values for tests or comparisons.
Edit active-noise parameters near the top of each file.
Current fixed-instantaneous-variance choice:

```matlab
eta_std = 0.5;
sigma_active = eta_std*sqrt(2/t_corr);
```
Fixed-integrated-strength alternative:

```matlab
noise_amplitude = 1.25;
sigma_active = sqrt(2*noise_amplitude)/t_corr;
```

Generated CSV files are ignored by Git.
