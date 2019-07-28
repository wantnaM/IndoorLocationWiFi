% clear; close all; clc;
% fc = 5.63e9;
% c = 3e8;
% deltaRange = [-c/2/fc, c/2/fc];
% M = 3;
% T = 2;
% d = 0.06;
% do_second_iter = 0;
% delayRange = [0 25e-9];
% SubCarrInd = -15:15;
% N = length(SubCarrInd);
% fgap = 312.5e3;
function [aTot,GridStart,GridSpacing, delayGridValue, u_sGridValue, deltaGridValue] = gridVecBackscatter(deltaRange, M, T, d, fc, c, do_second_iter, delayRange, SubCarrInd, N, fgap, GridPts, MaxAngle, MinAngle, generateAtot)
% do_second_iter = 0;

% Define angles at which MUSIC ?spectrum? will be computed
GridStart = [delayRange(1) MinAngle deltaRange(1)];
GridStop = [delayRange(2) MaxAngle deltaRange(2)];

% GridSpacing = (GridStop-GridStart)./(GridPts-ones(size(GridPts)));
GridSpacing = (GridStop-GridStart)./max(1, GridPts-ones(size(GridPts)));

angleStart = GridStart(2);
delayStart = GridStart(1);
deltaStart = GridStart(3);
angleDiff = GridSpacing(2);
delayDiff = GridSpacing(1);
deltaDiff = GridSpacing(3);
% siz = GridPts;
numGridPoints = prod(GridPts);

[delayIndices,angleIndices,deltaIndices] = ind2sub(GridPts,1:numGridPoints);
delayGridValue = delayStart + (delayIndices-1)*delayDiff;
angleGridValue = (angleStart + (angleIndices-1)*angleDiff)*pi/180;
u_sGridValue = (d*fc/c)*sin(angleGridValue);
deltaGridValue = deltaStart + (deltaIndices-1)*deltaDiff;

aTot = [];
% aTot = zeros(N*M*T,numGridPoints);
% for gridEntry = 1:numGridPoints
%     aTot(:,gridEntry) = gridSampleBackscatter(fc, T, deltaGridValue(gridEntry), M, u_sGridValue(gridEntry), c, SubCarrInd, fgap, delayGridValue(gridEntry) );
% end
if ~generateAtot
    load('aTotSaveSpotfi');
else
    aTot = zeros(N*M*T,numGridPoints);
    for gridEntry = 1:numGridPoints
        aTot(:,gridEntry) = gridSampleBackscatter(fc, T, deltaGridValue(gridEntry), M, u_sGridValue(gridEntry), c, SubCarrInd, fgap, delayGridValue(gridEntry) );
    end
    if generateAtot == 2
        save('aTotSaveSpotfi','aTot')
    end
end