%% Initialization (mostly the same for all examples here)
useGPU = true;

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

% Lattice parameters defination
Wavelength = 1064*10^-9;
CrossingAngle = 15/180*pi;
D =Wavelength/(2*sin(CrossingAngle/2)); % Lattice period 4 micons
Erecoil = hbar^2*(pi/D)^2/(2*m_Li6mol); % recoil energy 250 Hz.
LatticePotential = 50; % lattice depth in unit of recoil energy
V = @(X,Y,Z) 0.5*(X.^2 ); % function representing the external potential (must have 3 arguments)


xmax = 100*D/rscale; %  (50 micrometers)  grid boundaries [-xmax,xmax]
N=100000; % number of grid points for coordinate grid

if useGPU
    x = gpuArray.linspace(-xmax,xmax,N);
    grid = grid1d(x);
else
    grid = grid1d(xmax,N); % create a two-dimensional NxN grid object
end

task = GPEtask(grid,V); % initialize the GPE solver
task.Ntotal = 3000/2; % wave function normalization (total number of atoms)
% task.g = 4*pi*aLi/rscale/sqrt(2*pi/omegay*omegax); % nonlinear interaction constant
% task.g = 0;
task.g = 16*hbar^2* add_Li/(3*m_Li6mol* R_TF_y*R_TF_x); % nonlinear interaction constant
task.Xscale = rscale;
task.Tscale = tcoef;
task.LatticeConstant = D;

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
Deltat = 200*10^-6/tcoef; % total duration for simulation; (20 mico-seconds)
tstep = 1*10^-9/tcoef; % time step for calculation; (1 ns)
steps_int = 10;

% function representing the external potential (must have 3 arguments)
Vlattice = @(X,Y,Z) LatticePotential*Erecoil/(hbar*omegax)*cos((pi/D)*X*rscale).^2; 
task.UpdatePotential(Vlattice);
task.history.MomentumDensity = [];
task.user_callback = @FFTanalysis;
steps_ext = task.current_iter+ Deltat/tstep/steps_int;
task.solve_split(tstep,steps_int,steps_ext); % run the calculation

plot((1:steps_ext)*tstep*tcoef*10^6, task.history.MomentumDensity')
xlabel('Time/\mus')
ylabel("Nomalized Proportion")
title(['Lattice Depth: ',num2str(LatticePotential),'E_{recoil}'])
legend('-2hk','-hk','0','hk','2hk')
%% switch off the potential and TOF
% V = @(X,Y,Z) 0; 
% taskTOF = GPEtask(grid,V); % initialize the GPE solver
% taskTOF.Ntotal = 3000/2; % wave function normalization (total number of atoms)
% % task.g = 4*pi*aLi/rscale/sqrt(2*pi/omegay*omegax); % nonlinear interaction constant
% taskTOF.g = 0;
% taskTOF.Xscale = rscale;
% taskTOF.Tscale = tcoef;
% taskTOF.init_state = task.current_state;
% taskTOF.show_image = 1;
% 
% Deltat = 50000*10^-6/tcoef;
% tstep = 10000*10^-9/tcoef; % time step for calculation
% 
% 
% steps_ext = taskTOF.current_iter+ Deltat/tstep/steps_int;
% taskTOF.solve_split(tstep,steps_int,steps_ext); % run the calculation



%% fft analysis of the wavefunction
function ReturnText = FFTanalysis(task)

% plot the figure or not
PlotFigure = 0;
n = task.grid.nx;
NmaxP = 2;

Phi_FFT =abs( fftshift( fft(task.current_state))).^2/n;
if PlotFigure
    % Plot the FFT of phi
    hold off
    plot((-n/2:n/2-1),abs(Phi_FFT).^2)
    title('FFT analysis of Phi')
    xlim([-200 200])
    drawnow
end

hbarK =round( task.grid.weight*n / task.LatticeConstant * task.Xscale);


% derive the amplitude for the momentum
MomentumDesity = Phi_FFT((-NmaxP:NmaxP).*hbarK + ceil(n/2))  ;
MomentumDesity = MomentumDesity./sum(MomentumDesity);
task.history.MomentumDensity = [task.history.MomentumDensity;MomentumDesity];
ReturnText = ['0 mode Proportion: ',num2str(MomentumDesity(NmaxP+1),3) ];
end