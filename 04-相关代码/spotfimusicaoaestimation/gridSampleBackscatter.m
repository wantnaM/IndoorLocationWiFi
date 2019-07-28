function steeringVec = gridSampleBackscatter(fc, T, deltak, M, u_s, c, SubCarrInd, fgap, delay )
% INPUTS:
% fc: center frequency
% T = number of time instants with equal displacements
% deltak = dispacement along the thetak direction of the transmitter

aoaSteering = exp(-1i*2*pi*u_s*(-(M-1)/2:(M-1)/2)');
aoaSteeringInv = exp(-1i*2*pi*fc*(0:(T-1))'*deltak/c);
delaySteering = exp(-1i*2*pi*(SubCarrInd(:))*fgap*delay);

steeringVecDelayAoATmp = delaySteering*aoaSteering.';
steeringVecDelayAoA = steeringVecDelayAoATmp(:);
steeringVecTmp = steeringVecDelayAoA*aoaSteeringInv.';
steeringVec = steeringVecTmp(:);
