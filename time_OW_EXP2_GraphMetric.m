clear all
clc

maxNumCompThreads(1);

typeGG = 'RandLLE'; % log-linear #edges
% typeGG = 'RandSLE'; % sqrt-linear #edges

dsName = 'amazon';
maxKC = 1000;

load([dsName '_' num2str(maxKC) '_' typeGG '_Graph.mat']);

GM = zeros(nGG, nGG);
disp('compute the ground graph metric');
tic
for ii = 1:(nGG-1)
    [~, TRD_II, ~] = shortestpathtree(GG, ii, [(ii+1):nGG], 'OutputForm', 'cell');  
    GM(ii, (ii+1):nGG) = TRD_II; 
    GM((ii+1):nGG, ii) = TRD_II';
end
runTime_GroundGM = toc;
 
% histogram
XX_ID_vec = zeros(N, nGG);
tic
for ii = 1:N
    % WW{ii}
    tmpWW = WW{ii}/sum(WW{ii}); % normalization
    tmpXX_ID = XX_ID{ii};
    
    XX_ID_vec(ii, tmpXX_ID) = tmpWW'; 
end
runTime_Hist = toc;

% ================
nPair = 10000;
ff = load([dsName '_ID' num2str(nPair) '.mat']);
% ID: Nx2
ID = ff.ID;
DD_OW = zeros(nPair,1);

phi = @(X) exp(X.^2) - 1;
invphi = @(y) sqrt(log(y+1));
epsilon = 0.1;

tic
% compute the OT
for iiID = 1:nPair
    
    ii = ID(iiID, 1);
    jj = ID(iiID, 2);

    if mod(iiID, 5) == 0
        disp(['...' num2str(iiID)]);
    end

    % preprocessing
    tmpALL = XX_ID_vec(ii, :) + XX_ID_vec(jj, :); 
    idNZ = find(tmpALL > 0);
    
    tmpII = XX_ID_vec(ii, idNZ);
    tmpJJ = XX_ID_vec(jj, idNZ);
    GMIJ = GM(idNZ, idNZ);
        
    % ----------------
    DD_OW(iiID) = OrliczWasserstein(phi, invphi, tmpII', tmpJJ', GMIJ, epsilon);

end
runTime_Dist = toc;

runTime_Dist_ALL = runTime_Dist + runTime_GroundGM + runTime_Hist;

outName = [dsName '_Time_OW_EXP2_' num2str(maxKC) '_' typeGG '_' num2str(nPair) '.mat'];

avgRunTime = sum(runTime_Dist_ALL)/nPair;

save(outName, 'DD_OW', ...
     'runTime_Dist', 'runTime_GroundGM', 'runTime_Hist', 'runTime_Dist_ALL', ...
     'nPair', 'epsilon');
    
disp('FINISH !!!');


