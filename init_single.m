%Calculate the stationary state of GPE with Imaginary Time Propagation method
%
% Initialize basic parameters
%   Box size
L = single(40);
Lz = single(5);
%   Grid size
N=single(128);
Nz=single(128);
%   Physical parameters
g = single(0.0188);
NN=single(sqrt(6e5)); % square root of number of particles
niter = single(500); % number of ITP itrerations
dt = single(0.004); % time step

%   Initialize grid
r = linspace(-L/2,L/2,N);
rz = linspace(-Lz/2,Lz/2,Nz);
h = r(2)-r(1);
hz = rz(2)-rz(1);
k = [ (0:N/2)*2*pi/L -(N/2-1:-1:1)*2*pi/L];
kz = [ (0:Nz/2)*2*pi/Lz -(Nz/2-1:-1:1)*2*pi/Lz];


%% Initialize the potential V, initial condition phi0, momentum array kk and other necessary arrays

[X,Y,Z] = meshgrid(r,r,rz);
[KX,KY,KZ] = meshgrid(k,k,kz);
r0 = single(10.45);
om = single(0.5*4.88^2);
V = gpuArray(om*Z.^2 + 0.5*(sqrt(X.^2+Y.^2)-r0).^2);
phi0 = exp(-(om*Z.^2 + 0.5*(sqrt(X.^2+Y.^2)-r0).^2)/10);
phi = phi0*(NN/(sqrt(sum(sum(sum(abs(phi0).^2)))*h*h*hz)));
kk = gpuArray((KX.^2+KY.^2+KZ.^2)/2);
ekk = exp(-kk*dt);
%clear phi0 X Y Z KX KY KZ;
MU = zeros(niter,1,'gpuArray');
MU2 = zeros(niter,1,'gpuArray');
tmp = gpuArray(phi);

%% Do the main calculation
tic
tmp2 = real(tmp.*conj(tmp));
for i=1:niter
    tmp = exp(-(V + g*tmp2)*dt*0.5).*tmp;
    tmp = ifftn(ekk.*fftn(tmp));
    tmp = exp(-(V + g*tmp.*conj(tmp))*dt*0.5).*tmp;
	tmp2 = real(tmp.*conj(tmp));
    mu = NN/sqrt(sum(sum(sum(tmp2)))*h*h*hz);
    tmp=tmp*mu;
	tmp2 = tmp2*mu^2;
    MU(i) = mu;
    MU2(i) = sum(sum(sum(abs(conj(tmp).*ifftn(kk.*fftn(tmp)) + (V+g*tmp2).*tmp2))));
end
toc

%% Gather results and save
MU = 1/dt * log(gather(MU));
MU2 = gather(MU2)*h*h*hz./NN^2;
phi=gather(tmp);
save('MU','MU','MU2');
save('phi','phi');
imagesc(r,r,abs(phi(:,:,Nz/2)).^2);