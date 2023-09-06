%% fft analysis of the 2D wavefunction
% PlotFigure = 0;   plot the figure or not
% NmaxP = 100;      maximum mode of momentum to be analyzed
function ReturnText = FFTanalysis2D(task,PlotFigure,NmaxP,varargin)


nx = task.grid.nx;
ny = task.grid.ny;

if nargin == 4
    Nxfft = varargin{1};
    Nyfft = varargin{1};

else
    Nxfft = nx;
    Nyfft = ny;
end


if nargin == 5
    Xdirection = varargin{3};
    Ydirection = varargin{4};
else
    Xdirection = [1 0];
    Ydirection = [0 1];
end

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




% derive the amplitude for the momentum
SelectedMomentumX = (-NmaxP:NmaxP).*hbarKx*2 + ceil(Nxfft/2)+1;
SelectedMomentumY = (-NmaxP:NmaxP).*hbarKy*2 + ceil(Nyfft/2)+1;
MomentumDesity = Phi_FFT(SelectedMomentumX,SelectedMomentumY)  ;
% MomentumDesity = MomentumDesity./sum(MomentumDesity);
task.history.MomentumDensity{task.current_iter+1} = MomentumDesity;
ReturnText = ['0 mode Proportion: ',num2str(MomentumDesity(NmaxP+1,NmaxP+1)./sum(MomentumDesity(:)),3) ];
end
