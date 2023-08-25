%% Initialization (mostly the same for all examples here)

units % load some common physical constants
m_Li6mol = 2*m_Li6;
aLi_12 = 486*a0;
add_Li = aLi_12*0.6;
omegax = 16*2*pi; % trap frequencies 
omegay = 74*2*pi;
omegaz = 68*2*pi;

R_TF_x = 41*10^-6;
R_TF_y = 10*10^-6;

tcoef = 1/(omegax); % time scale
rscale = sqrt(hbar/m_Li6mol/omegax); % spatial scale

xmax = 25*D/rscale; %  (15 micrometers)  grid boundaries [-xmax,xmax]
N=25*500; % number of grid points for coordinate grid
grid = grid1d(xmax,N); % create a two-dimensional NxN grid object

% Lattice parameters defination
Wavelength = 1064*10^-9;
CrossingAngle = 15/180*pi;
D = Wavelength/(2*sin(CrossingAngle/2)); % Lattice period 4 micons
Erecoil = hbar^2*(pi/D)^2/(2*m_Li6mol); % recoil energy 250 Hz.
V = @(X,Y,Z) 0.5*(X.^2 ); % function representing the external potential (must have 3 arguments)

task = GPEtask(grid,V); % initialize the GPE solver
task.Ntotal = 3000/2; % wave function normalization (total number of atoms)
% task.g = 4*pi*aLi/rscale/sqrt(2*pi/omegay*omegax); % nonlinear interaction constant
% task.g = 0;
task.g = 16*hbar^2* add_Li/(3*m_Li6mol* R_TF_y*R_TF_x); % nonlinear interaction constant
task.Xscale = rscale;
task.Tscale = tcoef;
%% Stationary state calculation

tstep = 0.002; % initial time step for imaginary time evolution
acc = 1e-5; % desired accuracy
phi = task.groundstate_itp(tstep,acc); % phi will contain calculated stationary state

%% Initialize the harmonic potential

task.init_state = phi;
% tstep = 10*10^-9/tcoef; % time step for calculation
% steps_int = 50; % number of (internal) time steps in one (external) processing step
% steps_ext = 1; % number of (external) processing steps
task.show_image = 0; % show the solution on every precessing step
% task.solve_split(tstep,steps_int,steps_ext); % run the calculation


%% switch on lattice for Deltat
Deltat = 20*10^-6/tcoef;
tstep = 1*10^-9/tcoef; % time step for calculation
steps_int = 10;
% function representing the external potential (must have 3 arguments)
Vlattice = @(X,Y,Z) 500*Erecoil/(hbar*omegax)*cos((pi/D)*X*rscale).^2; 
task.UpdatePotential(Vlattice);
task.history.MomentumDensity = [];
task.user_callback = @FFTanalysis;
steps_ext = task.current_iter+ Deltat/tstep/steps_int;
task.solve_split(tstep,steps_int,steps_ext); % run the calculation


%% switch off the potential and TOF
V = @(X,Y,Z) 0; 
taskTOF = GPEtask(grid,V); % initialize the GPE solver
taskTOF.Ntotal = 3000/2; % wave function normalization (total number of atoms)
% task.g = 4*pi*aLi/rscale/sqrt(2*pi/omegay*omegax); % nonlinear interaction constant
taskTOF.g = 0;
taskTOF.Xscale = rscale;
taskTOF.Tscale = tcoef;
taskTOF.init_state = task.current_state;
taskTOF.show_image = 1;

Deltat = 50000*10^-6/tcoef;
tstep = 10000*10^-9/tcoef; % time step for calculation


steps_ext = taskTOF.current_iter+ Deltat/tstep/steps_int;
taskTOF.solve_split(tstep,steps_int,steps_ext); % run the calculation



%% fft analysis of the wavefunction
function MomentumDesity = FFTanalysis(task)

% Plot the FFT of phi
hold off
n = length(task.current_state);
Phi_FFT = fftshift( fft(task.current_state));
plot((-n/2:n/2-1),abs(Phi_FFT).^2)
title('FFT analysis of \phi')
xlim([-200 200])
drawnow

MomentumDesity = [];
% derive the 
end