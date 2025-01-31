function phi = solve_split(task,ddt,niter_inner,niter_outer)
% solve_split - Calculate the dynamics of GPE with the split-step method.
%
%  Usage :
%    phi = task.solve_split(ddt,niter_inner,niter_outer)
%  Input
%    ddt         :  evolution time step
%    niter_inner :  number of internal iterations between callbacks
%    niter_outer :  number of external iterations
%  Output
%    phi       :  final state

task.dispstat('','init');
tic;
grid = task.grid;
% VV = task.getVtotal(0);
g = task.g;
omega = task.omega;
task.history.tstep = ddt;
task.history.steps_int = niter_inner;
task.history.steps_ext = niter_outer;

if(task.Ntotal > 0)
    NN0 = task.Ntotal;
else
    NN0 = grid.integrate(abs(task.init_state).^2);
end
NNN=NN0;
gam = task.gamma;
n_cn = task.n_crank;
n_rec = task.n_recalc;
tau = task.decay_rate;
start = task.current_iter;

dt = ddt*1i/(1+1i*gam);
ekk = exp(-grid.kk*dt);

if(start>0)
    phi = task.current_state;
else
    phi = task.init_state;
    task.ext_callback(phi,0,(ddt*n_rec),task.mu_init,task.Ntotal);
end

tmp2 = real(phi.*conj(phi));
muc = real(grid.inner(phi,task.applyham(phi)))./NN0;
if(task.mu_init > 0)
    mu = task.mu_init;
else
    mu = muc;
end
dt_outer = ddt*niter_inner;
% main BIG cycle starts here
for j=start+1:niter_outer
    time=(j-1)*dt_outer;
    for jj=1:niter_inner/n_rec
        
        time2=time+(jj-1)*ddt*n_rec;
        VV = task.getVtotal(time2)-mu;
        phi = exp(( - VV - g.*tmp2)*dt*0.5).*phi;
        % main SMALL cycle starts here
        for i=1:n_rec
            phi = grid.ifft(ekk.*grid.fft(phi));
            if(omega ~= 0)
                lphi = phi;
                for ii = 1:n_cn
                    lphi = phi + dt*omega*grid.lz(lphi);
                    lphi = 0.5*(phi+lphi);
                end
                phi = phi + dt*omega*grid.lz(lphi);
            end
            
            phi = exp(( - VV - g.*phi.*conj(phi))*dt).*phi;
        end
        
        phi = exp((VV + g.*phi.*conj(phi))*dt*0.5).*phi;
        tmp2 = real(phi.*conj(phi));
        
        ncur = grid.integrate(tmp2);
        if(gam>0 && task.mu_init == 0)
            if(tau >0)
                NNN = NN0*exp(-time2/tau);
            end                
            phi = phi*sqrt(NNN/ncur);
            tmp2 = tmp2*(NNN/ncur);
            muc = real(grid.inner(phi,task.applyham(phi,time2)))/NNN;
            mu = muc;
        else
            muc = real(grid.inner(phi,task.applyham(phi,time2)))/ncur;
        end
    end
    task.ext_callback(phi,j,(time2+ddt*n_rec),muc,ncur);
    
end

end
