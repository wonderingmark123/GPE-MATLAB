%% Initialization (mostly the same for all examples here)
% The following is the simulation code for 
% Liang, Qi, et al. "Diffraction of strongly interacting molecular Bose-Einstein condensate from standing wave light pulses." 
% SciPost Physics 12.5 (2022): 154.



units % load some common physical constants
omegar = 150*2*pi; % trap frequencies (150Hz and 600Hz)
omegaz = 600*2*pi;
tcoef = 1/(omegar); % time scale
rscale = sqrt(hbar/mRB/omegar); % spatial scale

xmax = 15e-6/rscale; %  (15 micrometers)  grid boundaries [-xmax,xmax]
N=256; % number of grid points for coordinate grid
grid = grid2d(xmax,N,xmax,N); % create a two-dimensional NxN grid object

V = @(X,Y,Z) 0.5*(X.^2 + Y.^2); % function representing the external potential (must have 3 arguments)
task = GPEtask(grid,V); % initialize the GPE solver
task.Ntotal = 1.0e5; % wave function normalization (total number of atoms)
task.g = 4*pi*aRB/rscale/sqrt(2*pi/omegaz*omegar); % nonlinear interaction constant

%% Stationary state calculation

tstep = 0.02; % initial time step for imaginary time evolution
acc = 1e-6; % desired accuracy
phi = task.groundstate_itp(tstep,acc); % phi will contain calculated stationary state

%% Initialize the harmonic potential

task.init_state = phi;
tstep = 0.0004; % time step for calculation
steps_int = 50; % number of (internal) time steps in one (external) processing step
steps_ext = 10; % number of (external) processing steps
task.show_image = 1; % show the solution on every precessing step
task.solve_split(tstep,steps_int,steps_ext); % run the calculation


%% switch on lattice for Deltat
Deltat = 8;
k = 1;

Vlattice = @(X,Y,Z) 100*cos(k*X); % function representing the external potential (must have 3 arguments)
task.UpdatePotential(Vlattice);

steps_ext = task.current_iter+ Deltat/tstep/steps_int;
task.solve_split(tstep,steps_int,steps_ext); % run the calculation

