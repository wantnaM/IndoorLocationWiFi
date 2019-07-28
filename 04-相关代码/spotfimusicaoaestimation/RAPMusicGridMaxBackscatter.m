function [delayFromMusic, angleFromMusic, deltaFromMusic, maxCorr, music_spectrum_plot] = RAPMusicGridMaxBackscatter(aTot,GridStart,GridSpacing,GridPts,Qn,Qs,fc,fgap,d,K,L,delayFromMusicPrev,angleFromMusicPrev, deltaFromMusicPrev, ...
                                                    SubCarrInd, deltaGridValue, u_sGridValue, delayGridValue, T, c, do_second_iter, seconditerGridPts, useNoise)
maxCorr = [];
if ~useNoise
    if ~isempty(delayFromMusicPrev)
        Ahat = [];
        for compNo = 1:length(delayFromMusicPrev)
            u_s = (d*fc/c)*sind(angleFromMusicPrev(compNo));
            Ahat(:,compNo) = gridSampleBackscatter(fc, T, deltaFromMusicPrev(compNo), K, u_s, c, SubCarrInd(1:L), fgap, delayFromMusicPrev(compNo) );
        end
        PerpAhat = eye(size(Ahat,1)) - Ahat*pinv(Ahat'*Ahat)*Ahat';
    else
        PerpAhat = eye(size(Qs,1),'like',2+1i);
    end
    Qs = PerpAhat*Qs;
    Ps = (Qs*Qs');

    numGridPts = prod(GridPts);

    % achieving the music spectrum without for loop by only matrix
    % computations
    delayConsider = GridStart(1) + (0:GridPts(1)-1)*GridSpacing(1);
    u_sConsider = (d*fc/c)*sind( GridStart(2) + (0:GridPts(2)-1)*GridSpacing(2) );
    deltaConsider = GridStart(3) + (0:GridPts(3)-1)*GridSpacing(3);
    
    delaySteeringMat = exp(-1i*2*pi*(SubCarrInd(1:L).')*fgap*delayConsider);
    aoaSteeringMat = exp(-1i*2*pi*((-(K-1)/2:(K-1)/2)')*u_sConsider);
    aoaSteeringInvMat = exp(-1i*2*pi*(fc/c)*((0:(T-1))')*deltaConsider);

    thetaTauMat = kron(aoaSteeringMat, delaySteeringMat);
    thetaTauDeltaMat = kron(aoaSteeringInvMat, thetaTauMat);

    % we have now created steering vectors where each column is a steering
    % vector.
    PerpAhat_a = PerpAhat*thetaTauDeltaMat; % modifying each steering vector. These are the a1's.
    PsA = Ps*PerpAhat_a;
    music_spec_num = sum((PerpAhat_a').*(PsA.'),2);            % now conjugate of each a1 should be dot producted with columns of PsA
    music_spec_den = sum((PerpAhat_a').*(PerpAhat_a.'),2); % now conj of each a1 should be dot producted with a1
    music_spectrum = abs(music_spec_num./music_spec_den);
    
%     % previous way of computing MUSIC spectrum by for looping over all
%     % possible grid values
%     music_spectrum = zeros(numGridPts,1);
%     % Use PARFOR
%     for loopVar = 1:numGridPts
%         if ~isempty(aTot)
%             aTotTmp = aTot(:,loopVar); % use this if you have precomputed the grid
%         else
%             aTotTmp = gridSampleBackscatter(fc, T, deltaGridValue(loopVar), K, u_sGridValue(loopVar), c, SubCarrInd(1:L), fgap, delayGridValue(loopVar) );
%         end
%         a1 = PerpAhat*aTotTmp;
%         music_spectrum(loopVar) = (a1'*Ps*a1)/(a1'*a1); 
%     end
%     music_spectrum = abs(music_spectrum);
    
else
    % delay forms first dimension, angle second and displacement third
    delayConsider = GridStart(1) + (0:GridPts(1)-1)*GridSpacing(1);
    u_sConsider = (d*fc/c)*sind( GridStart(2) + (0:GridPts(2)-1)*GridSpacing(2) );
    deltaConsider = GridStart(3) + (0:GridPts(3)-1)*GridSpacing(3);
%     delayConsider = unique(delayGridValue, 'stable');
%     u_sConsider = unique(u_sGridValue, 'stable');
%     deltaConsider = unique(deltaGridValue, 'stable');
    
    delaySteeringMat = exp(-1i*2*pi*(SubCarrInd(1:L).')*fgap*delayConsider);
    aoaSteeringMat = exp(-1i*2*pi*((-(K-1)/2:(K-1)/2)')*u_sConsider);
    aoaSteeringInvMat = exp(-1i*2*pi*(fc/c)*((0:(T-1))')*deltaConsider);
    
    thetaTauMat = kron(aoaSteeringMat, delaySteeringMat);
    thetaTauDeltaMat = kron(aoaSteeringInvMat, thetaTauMat);
    
    PnA = (Qn*Qn')*thetaTauDeltaMat;
    thetaTauDeltaMatTrans = thetaTauDeltaMat';
    % we want only the diagonal terms of thetaTauDeltaMatTrans*PnA
    music_spectrum = sum(thetaTauDeltaMatTrans.*(PnA.'),2);

    music_spectrum = 1./abs(music_spectrum);
%     figure; imagesc(reshape(squeeze(music_spectrum), GridPts([1 3])))
%     figure; imagesc(reshape((music_spectrum), GridPts([1 2])))
end

%% Finding maximum of Music specturm and the corresponding Grid values
if ~useNoise
    [maxCorr,loopVarMax] = max(music_spectrum); % returns the first occurence of the maximum
    % if doGPU
    %     maxCorr = gather(maxCorr);
    %     loopVarMax = gather(loopVarMax);
    % end
    [delay_idx,angle_idx, delta_idx] = ind2sub(GridPts,loopVarMax); 
    delayFromMusic = GridStart(1) + (delay_idx-1)*GridSpacing(1);
    angleFromMusic = GridStart(2) + (angle_idx-1)*GridSpacing(2);
    deltaFromMusic = GridStart(3) + (delta_idx-1)*GridSpacing(3);
else
    music_spectrum = reshape(music_spectrum,GridPts(1),GridPts(2),GridPts(3));
    BW = imregionalmax(music_spectrum);
    [delay_idx, angle_idx, delta_idx] = ind2sub(size(BW), find(BW));
    topPeakIndices = topk(music_spectrum(BW), min(size(Qs,2), length(find(BW(:)))));
    delayFromMusic = GridStart(1) + (delay_idx(topPeakIndices)-1)*GridSpacing(1);
    angleFromMusic = GridStart(2) + (angle_idx(topPeakIndices)-1)*GridSpacing(2);
    deltaFromMusic = GridStart(3) + (delta_idx(topPeakIndices)-1)*GridSpacing(3);
    maxCorr = max(music_spectrum);
end

music_spectrum = reshape(music_spectrum,GridPts(1),GridPts(2),GridPts(3));
music_spectrum_plot = music_spectrum;

%% for second iteration of MUSIC

if do_second_iter  
    
    delayFirstIter = delayFromMusic;
    angleFirstIter = angleFromMusic;
    deltaFirstIter = deltaFromMusic;
    
    for iComp = 1:length(delayFirstIter)
        GridPts = seconditerGridPts;
        delayRange = [delayFirstIter(iComp)-GridSpacing(1) delayFirstIter(iComp)+GridSpacing(1)];
        angleRange = [angleFirstIter(iComp)-GridSpacing(2) angleFirstIter(iComp)+GridSpacing(2)];
        deltaRange = [deltaFirstIter(iComp)-GridSpacing(3) deltaFirstIter(iComp)+GridSpacing(3)];
        [GridStart, GridSpacing, delayGridValue, u_sGridValue, deltaGridValue] = gridParamsBackscatter(GridPts, angleRange, deltaRange, d, fc, c, delayRange);

        if ~useNoise
            numGridPts = prod(GridPts);

            % achieving the music spectrum without for loop by only matrix
            % computations
            delayConsider = GridStart(1) + (0:GridPts(1)-1)*GridSpacing(1);
            u_sConsider = (d*fc/c)*sind( GridStart(2) + (0:GridPts(2)-1)*GridSpacing(2) );
            deltaConsider = GridStart(3) + (0:GridPts(3)-1)*GridSpacing(3);

            delaySteeringMat = exp(-1i*2*pi*(SubCarrInd(1:L).')*fgap*delayConsider);
            aoaSteeringMat = exp(-1i*2*pi*((-(K-1)/2:(K-1)/2)')*u_sConsider);
            aoaSteeringInvMat = exp(-1i*2*pi*(fc/c)*((0:(T-1))')*deltaConsider);

            thetaTauMat = kron(aoaSteeringMat, delaySteeringMat);
            thetaTauDeltaMat = kron(aoaSteeringInvMat, thetaTauMat);

            % we have now created steering vectors where each column is a steering
            % vector.
            PerpAhat_a = PerpAhat*thetaTauDeltaMat; % modifying each steering vector. These are the a1's.
            PsA = Ps*PerpAhat_a;
            music_spec_num = sum((PerpAhat_a').*(PsA.'),2);            % now conjugate of each a1 should be dot producted with columns of PsA
            music_spec_den = sum((PerpAhat_a').*(PerpAhat_a.'),2); % now conj of each a1 should be dot producted with a1
            music_spectrum = abs(music_spec_num./music_spec_den);

            % % %         Previous way of doing it through far loop
            % music_spectrum = zeros(numGridPts,1);
            % for loopVar = 1:numGridPts
            %     aTotTmp = gridSampleBackscatter(fc, T, deltaGridValue(loopVar), K, u_sGridValue(loopVar), c, SubCarrInd(1:L), fgap, delayGridValue(loopVar) );
            %     a1 = PerpAhat*aTotTmp;
            %     music_spectrum(loopVar) = (a1'*Ps*a1)/(a1'*a1); 
            % end
            % music_spectrum = abs(music_spectrum);
        else
            % delay forms first dimension, angle second and displacement third
            delayConsider = GridStart(1) + (0:GridPts(1)-1)*GridSpacing(1);
            u_sConsider = (d*fc/c)*sind( GridStart(2) + (0:GridPts(2)-1)*GridSpacing(2) );
            deltaConsider = GridStart(3) + (0:GridPts(3)-1)*GridSpacing(3);
            % delayConsider = unique(delayGridValue, 'stable');
            % u_sConsider = unique(u_sGridValue, 'stable');
            % deltaConsider = unique(deltaGridValue, 'stable');

            delaySteeringMat = exp(-1i*2*pi*(SubCarrInd(1:L).')*fgap*delayConsider);
            aoaSteeringMat = exp(-1i*2*pi*((-(K-1)/2:(K-1)/2)')*u_sConsider);
            aoaSteeringInvMat = exp(-1i*2*pi*(fc/c)*((0:(T-1))')*deltaConsider);

            thetaTauMat = kron(aoaSteeringMat, delaySteeringMat);
            thetaTauDeltaMat = kron(aoaSteeringInvMat, thetaTauMat);

            PnA = (Qn*Qn')*thetaTauDeltaMat;
            thetaTauDeltaMatTrans = thetaTauDeltaMat';
            % we want only the diagonal terms of thetaTauDeltaMatTrans*PnA
            music_spectrum = sum(thetaTauDeltaMatTrans.*(PnA.'),2);

            music_spectrum = 1./abs(music_spectrum); 
        end

        if ~useNoise
            [maxCorrTmp,loopVarMaxTmp] = max(music_spectrum); 
            maxCorr = maxCorrTmp(1);
            loopVarMax = loopVarMaxTmp(1);
            [delay_idx,angle_idx, delta_idx] = ind2sub(GridPts,loopVarMax); 
            delayFromMusic = GridStart(1) + (delay_idx-1)*GridSpacing(1);
            angleFromMusic = GridStart(2) + (angle_idx-1)*GridSpacing(2);
            deltaFromMusic = GridStart(3) + (delta_idx-1)*GridSpacing(3);
        else
            maxCorr = max(music_spectrum);
            music_spectrum = reshape(music_spectrum,GridPts(1),GridPts(2),GridPts(3));
            BW = imregionalmax(music_spectrum);
            [delay_idx, angle_idx, delta_idx] = ind2sub(size(BW), find(BW));
            topPeakIndices = topk(music_spectrum(BW), 1);
            delayFromMusic(iComp) = GridStart(1) + (delay_idx(topPeakIndices)-1)*GridSpacing(1);
            angleFromMusic(iComp) = GridStart(2) + (angle_idx(topPeakIndices)-1)*GridSpacing(2);
            deltaFromMusic(iComp) = GridStart(3) + (delta_idx(topPeakIndices)-1)*GridSpacing(3);
        end
    end
end

