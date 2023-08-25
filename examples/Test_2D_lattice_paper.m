%% Initialization (mostly the same for all examples here)

units % load some common physical constants
m_Li6mol = 2*m_Li6;
aLi = 486*a0;
omegax = 16*2*pi; % trap frequencies 
omegay = 74*2*pi;
omegaz = 68*2*pi;

tcoef = 1/(omegax); % time scale
rscale = sqrt(hbar/mRB/omegax); % spatial scale

xmax = 25*D/rscale; %  (15 micrometers)  grid boundaries [-xmax,xmax]
N=256; % number of grid points for coordinate grid
grid = grid2d(xmax,N,xmax,N); % create a two-dimensional NxN grid object

% Lattice parameters defination
Wavelength = 1064*10^-9;
CrossingAngle = 15/180*pi;
D = Wavelength/(2*sin(CrossingAngle/2)); % Lattice period 4 micons
Erecoil = hbar^2*(pi/D)^2/(2*m_Li6mol); % recoil energy 250 Hz.
V = @(X,Y,Z) 0.5*(X.^2 + Y.^2); % function representing the external potential (must have 3 arguments)

task = GPEtask(grid,V); % initialize the GPE solver
task.Ntotal = 3000/2; % wave function normalization (total number of atoms)
task.g = 4*pi*aLi/rscale/sqrt(2*pi/omegay*omegax); % nonlinear interaction constant
% task.g = 0;
%% Stationary state calculation

tstep = 0.02; % initial time step for imaginary time evolution
acc = 1e-5; % desired accuracy
phi = task.groundstate_itp(tstep,acc); % phi will contain calculated stationary state

%% Initialize the harmonic potential

task.init_state = phi;
tstep = 100*10^-9/tcoef; % time step for calculation
steps_int = 50; % number of (internal) time steps in one (external) processing step
steps_ext = 10; % number of (external) processing steps
task.show_image = 1; % show the solution on every precessing step
task.solve_split(tstep,steps_int,steps_ext); % run the calculation


%% switch on lattice for Deltat
Deltat = 50000*10^-6/tcoef;


% function representing the external potential (must have 3 arguments)
Vlattice = @(X,Y,Z) 500*Erecoil/(hbar*omegax)*cos((pi/D)*X*rscale); 
task.UpdatePotential(Vlattice);

steps_ext = task.current_iter+ Deltat/tstep/steps_int;
task.solve_split(tstep,steps_int,steps_ext); % run the calculation

%% switch off the potential and TOF
Deltat = 10000*10^-6/tcoef;


% function representing the external potential (must have 3 arguments)
Vlattice = @(X,Y,Z) 0; 
task.UpdatePotential(Vlattice);

steps_ext = task.current_iter+ Deltat/tstep/steps_int;
task.solve_split(tstep,steps_int,steps_ext); % run the calculation
