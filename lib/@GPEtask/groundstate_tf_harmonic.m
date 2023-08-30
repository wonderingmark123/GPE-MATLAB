function [phi,mu,R] = groundstate_tf_harmonic(task,a,wx,wy,wz,N,V)
% groundstate_tf_harmonic - Calculate the Thomas-Fermi stationary state in harmonic trap.
%
%  Usage :
%    [phi, mu] = task.groundstate_tf(eps)
%    [phi, mu] = task.groundstate_tf(eps,N,V,mul,mur)
%  Input
%    eps     :  desired accuracy (applied to chemical potential)
%    wx      :  Trap frequency in x direction unit is radial
%    wy      :  Trap frequency in y direction
%    wz      :  Trap frequency in z direction
%    V       :  trap potential 
%    N       :  particle number
%    a       :  interaction length

%  Output
%    phi      :  calculated stationary state
%    mu       :  chemical potential

% 

if(nargin < 7)
    if(nargin < 6)
        N = task.Ntotal;
    end
    V = task.getVtotal(0);
end




mu = 1/2*(15*N*a/task.Xscale).^(2/5) * (wx*wy*wz).^(1/3)/wx;
Rx = sqrt(2*mu);
Ry = sqrt(2*mu)*wx/wy;
Rz = sqrt(2*mu)*wx/wz;


vvv = (mu-V).*(V<mu);
phi = sqrt(vvv);


tmp2 = real(phi.*conj(phi));        
ncur = task.grid.integrate(tmp2);
phi = phi*sqrt(task.Ntotal/ncur);

R = [Rx,Ry,Rz];
end