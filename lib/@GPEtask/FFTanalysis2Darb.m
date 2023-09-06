%% fft analysis of the 2D wavefunction 
% get the results of arbitrary mementum.
% PlotFigure = 0;   plot the figure or not
% Px = [-1 0 1];      maximum mode of momentum to be analyzed
function ReturnText = FFTanalysis2Darb(task,PlotFigure,Px,Py)


nx = task.grid.nx;
ny = task.grid.ny;

Nxfft = nx;
Nyfft = ny;

hbarKx =round(abs(task.grid.x(1)) / task.LatticeConstant * task.Xscale*Nxfft/nx);
hbarKy =round(abs(task.grid.y(1)) / task.LatticeConstant * task.Xscale*Nyfft/ny);

% Nxfft = hbarKx * NmaxP *2 * scale;
% Nyfft = hbarKy * NmaxP *2 * scale;

Phi_FFT =abs( fftshift( fft2(task.current_state',Nxfft,Nyfft))).^2;
if PlotFigure
    % Plot the FFT of phi
    hold off
    %     subplot 211
    imagesc((Phi_FFT).^(0.1))
    title('FFT analysis of \Phi')
    colormap pink
    %     subplot 212
    %     plot(Phi_FFT(501,:))
    drawnow
end

maxPx = length(Px);
MomentumDesity = zeros(maxPx,1);

for i = 1:length(Px)

% derive the amplitude for the momentum
SelectedMomentumX =round( Px(i).*hbarKx*2 + ceil(Nxfft/2)+1);
SelectedMomentumY =round( Py(i).*hbarKy*2 + ceil(Nyfft/2)+1);
MomentumDesity(i) = Phi_FFT(SelectedMomentumX,SelectedMomentumY)  ;
% MomentumDesity = MomentumDesity./sum(MomentumDesity);

end
task.history.MomentumDensity{task.current_iter+1} = MomentumDesity;
ReturnText = '';
end