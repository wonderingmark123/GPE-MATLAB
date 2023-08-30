%% fft analysis of the wavefunction
% PlotFigure = 0;   plot the figure or not
% NmaxP = 100;      maximum mode of momentum to be analyzed
function ReturnText = FFTanalysis(task,PlotFigure,NmaxP)


n = task.grid.nx;


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
% MomentumDesity = MomentumDesity./sum(MomentumDesity);
task.history.MomentumDensity = [task.history.MomentumDensity;MomentumDesity];
ReturnText = ['0 mode Proportion: ',num2str(MomentumDesity(NmaxP+1)./sum(MomentumDesity),3) ];
end
