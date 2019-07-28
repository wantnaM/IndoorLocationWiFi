
function [Pn,Ps,Qn,Qs,EigenInfo] = GetQnBackscatter(X,EigDiffCutoff, nComps)

% % by singular value
% [U,D,~] = svd(X);
% % by eigenvalue
[Utmp,D] = eig(X*X');
D = abs(D);
[Dtmp,I] = sort(diag(D), 'descend');
D = diag(Dtmp);
U = Utmp(:,I);

% % % density based clustering:
% [class,~] = DBSCAN(diag(D),2);
% SignalEndIdx = find(class==-1, 1, 'last');

minMP = 2;
useMDL = 0;
useDiffMaxVal = 0; % Default value is 1. If set to 1, it considers only those eignevalues who are above a certain threshold when compared to the maximum eigenvalue

% % % MDL criterion based
MDL = [];
lambdaTot = diag(D);
subarraySize = size(X,1);
nSegments = size(X,2);
maxMultipath = length(lambdaTot); % previously used 6 for maximum number of multipath
for k = 1:maxMultipath
    MDL(k) = -nSegments*(subarraySize-(k-1))*log(geo_mean(lambdaTot(k:end))/mean(lambdaTot(k:end))) + 0.5*(k-1)*(2*subarraySize-k+1)*log(nSegments);
end
% % Another attempt to take the number of multipath as minimum of MDL
[~, SignalEndIdxTmp] = min(MDL);

% % TEMP changing so that small changes in trough of MDL values dont affect
% % the signalEndIdx value
% [~, max_mdl] = min(MDL);
% mdlConsider = MDL(1:max_mdl);
% [class_mdl, type_mdl] = DBSCAN(mdlConsider(:), 1, 0.05*(max(mdlConsider)-min(mdlConsider))); % used 0.05 before
% SignalEndIdxTmp = max_mdl; % consider the smallest MDL to be 1+number of components
% if type_mdl(end)~=-1    % keep the max MDL the same as previous if the end point is not part of any cluster
%     % if ~isempty(find(type_mdl~=-1, 1))  % keep MDL same if no point belongs to cluster, this is redudndant
%         SignalEndIdxTmp = find(class_mdl == class_mdl(end), 1);
%     % end
% end

if useMDL
    SignalEndIdx = max(SignalEndIdxTmp-1, 1);
else
    % % % Older way of finding the SignalEndIdx based on thresholding and
    % eigenvalue difference cutoff
    VecForSignalEnd = wkeep(diag(D),5,'l'); % latest last used is 5 % Previously used 8 which means allow upto 8 eigenvalues
    diag(D(1:floor(length(D)/2)));
    % SignalEndIdx = min(3,find(diff(db(VecForSignalEnd))<-EigDiffCutoff,1,'last'));
    % SignalEndIdx = length(find(VecForSignalEnd/VecForSignalEnd(1)>0.065));
    Criterion1 = diff(db(VecForSignalEnd))<=max(-EigDiffCutoff,min(diff(db(VecForSignalEnd))));
    Criterion3 = (VecForSignalEnd(1:end-1)/VecForSignalEnd(1)>0.03); %  Previously used 0.165 and right now using 0.065
    % % % previously used criterion
    % SignalEndIdx = find(Criterion1,1,'last');
    SignalEndIdx = find(Criterion1 & Criterion3,1,'last');
    if isempty(SignalEndIdx)
        SignalEndIdx = find(Criterion3,1,'last');
    end
end

SignalEndIdx = max(SignalEndIdx, minMP);
if ~isempty(nComps)
    SignalEndIdx = nComps; 
end

% % % % % % Plotting the figures to demnstrate difference in eigenvalues.
% figure(1); plot(diff(db(diag(D))),'d-')
% % % SignalEndIdx = 10;
% figure(2); plot((diag(D)),'d-')
% figure(3); plot(zeros(length(diag(D)),1), diag(D),'d');
% hold on;
% % plot([-1 1], 0.01*[1 1], 'r--')
% plot([-1 1], 0.1*D(1)*[1 1], 'r--')
% hold off;
% figure(4); plot(MDL, 'd-');
% pause(0.1)
% % sprintf('min criterion is %d', find(Criterion1 & Criterion3,1,'last'))
% sprintf('Noise space dimesion is %d',size(U,2)-SignalEndIdx)
% sprintf('Signal space dimesion is %d',SignalEndIdx)
% sprintf('min MDL index is %d', SignalEndIdxTmp)

Qn = U(:,SignalEndIdx+1:end);
Pn = (Qn*Qn');

Qs = U(:,1:SignalEndIdx);
Ps = (Qs*Qs');

EigenInfo = struct;
EigenInfo.Umatrix = U;
EigenInfo.nMPs = SignalEndIdx;
EigenInfo.singularValuesdB = db(diag(D));
% EigenInfo.CSI = ysenarr;
% varargout{1} = EigenInfo;
