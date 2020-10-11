%%
%     COURSE: Advanced neuroscience techniques
%    LECTURE: Spatial source separation
% Instructor: mikexcohen.com
%
%%


%%
%
% Three ways to make a covariance matrix
%
%

% load the dataset. dimensions are channels X time X repetitions (trials)
load multivar_data.mat

% we'll work with the trial-average here
data = mean(data,3);

datadim = size(data);

%% the method with loops

% initialize
covmat1 = zeros(datadim(1));

% double-loop over channels and compute dot product scaled by N-1
for chani=1:datadim(1)
    for chanj=1:datadim(1)
        
        % mean-center data
        subi = data(chani,:)
        subj = data(chanj,:)
        
        % compute covariance as dot product divided by N-1
        covmat1(chani,chanj) = 
    end % end chanj
end % end chani


%% using matrix multiplication

% first mean-center (over time!)
dataM = bsxfun(@minus,data,mean(data,2));

% matrix times its transpose, divided by N-1
covmat2 = dataM*dataM' / ;


%% using MATLAB's cov function

covmat3 = cov(data);

% always check the size to make sure it's the correct orientation:
% size(covmat3)

%% show all three in a plot

figure(1), clf

titles = {'Using loops';'Matrix multiplication';'MATLAB cov function'};


for i=1:3
    subplot(1,3,i)
    
    % use eval for dynamic variable naming
    eval([ 'imagesc(covmat' num2str(i) ')' ])
    axis square
    title(titles{i})
    xlabel('Channel'), ylabel('Channel')
end

%%% QUESTION: Do the three covariance matrices look the same?

%%

% load mat file containing EEG, leadfield and channel locations
load emptyEEG

% pick a dipole location in the brain
diploc = 109;


% normalize dipoles
lf.GainN = bsxfun(@times,squeeze(lf.Gain(:,1,:)),lf.GridOrient(:,1)') + bsxfun(@times,squeeze(lf.Gain(:,2,:)),lf.GridOrient(:,2)') + bsxfun(@times,squeeze(lf.Gain(:,3,:)),lf.GridOrient(:,3)');


% plot brain dipoles
figure(1), clf, subplot(221)
plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), 'o')
hold on
plot3(lf.GridLoc(diploc,1), lf.GridLoc(diploc,2), lf.GridLoc(diploc,3), 's','markerfacecolor','w','markersize',10)
rotate3d on, axis square, axis off
title('Brain dipole locations')


% Each dipole can be projected onto the scalp using the forward model. 
% The code below shows this projection from one dipole.
subplot(222)
topoplotIndie(lf.GainN(:,diploc), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
title('Signal dipole projection')


% Now we generate random data in brain dipoles.
% create 1000 time points of random data in brain dipoles
% (note: the '1' before randn controls the amount of noise)
dipole_data = 1*randn(length(lf.Gain),1000);

% add signal to second half of dataset
dipole_data(diploc,501:end) = 15*sin(2*pi*10*(0:499)/EEG.srate);

% project dipole data to scalp electrodes
EEG.data = lf.GainN*dipole_data;

% meaningless time series
EEG.times = (0:size(EEG.data,2)-1)/EEG.srate;

% plot the data from one channel
subplot(212), hold on, plot(.5,0,'HandleVisibility','off');
plot(EEG.times,dipole_data(diploc,:)/norm(dipole_data(diploc,:)),'linew',4)
plot(EEG.times,EEG.data(31,:)/norm(EEG.data(31,:)),'linew',2)
xlabel('Time (s)'), ylabel('Amplitude (norm.)')
legend({'Dipole';'Electrode'})

%% Create covariance matrices

% compute covariance matrix S (signal data)
tmpd = 
tmpd = 
covS = 

% compute covariance matrix R (reference data)
tmpd = % select data
tmpd = % mean-center
covR = % create covariance matrix


%%% plot the two covariance matrices
figure(2), clf

% S matrix
subplot(131)
imagesc(covS)
title('S matrix')
axis square, set(gca,'clim',[-1 1]*1e6)

% R matrix
subplot(132)
imagesc(covR)
title('R matrix')
axis square, set(gca,'clim',[-1 1]*1e6)

% R^{-1}S
subplot(133)
imagesc(inv(covR)*covS)
title('R^-^1S matrix')
axis square, set(gca,'clim',[-10 10])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                 %
%  Dimension compression via PCA  %
%                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PCA
[evecs,evals] = eig(); % what to take the eigendecomposition of?

% sort eigenvalues/vectors
[evals,sidx] = sort(diag(evals),'descend');
% DON'T FORGET: you need to sort the eigenvectors according to the
% eigenvalues sorting.



% plot the eigenspectrum
figure(3), clf
subplot(231)
plot(evals./max(evals),'s-','markersize',15,'markerfacecolor','k')
axis square
set(gca,'xlim',[0 20.5])
title('PCA eigenvalues')
xlabel('Component number'), ylabel('Power ratio (norm-\lambda)')


% component time series is eigenvector as spatial filter for data
comp_ts = evecs(:,1)'*EEG.data;


% normalize time series (for visualization)
dipl_ts = dipole_data(diploc,:) / norm(dipole_data(diploc,:));
comp_ts = comp_ts / norm(comp_ts);
chan_ts = EEG.data(31,:) / norm(EEG.data(31,:));


% plot the time series
subplot(212), hold on
plot(EEG.times,.3+dipl_ts,'linew',2)
plot(EEG.times,.15+chan_ts)
plot(EEG.times,comp_ts)
legend({'Truth';'EEG channel';'PCA time series'})
set(gca,'ytick',[])
xlabel('Time (a.u.)')


%% spatial filter forward model

% The filter forward model is what the source "sees" when it looks through the
% electrodes. It is obtained by passing the covariance matrix through the filter.
filt_topo = evecs(:,1);

% Eigenvector sign uncertainty can cause a sign-flip, which is corrected for by 
% forcing the largest-magnitude projection electrode to be positive.
[~,se] = max(abs( filt_topo ));
filt_topo = filt_topo * sign(filt_topo(se));


% plot the maps
subplot(232)
topoplotIndie(lf.GainN(:,diploc), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
title('Truth topomap')

subplot(233)
topoplotIndie(filt_topo,EEG.chanlocs,'electrodes','numbers','numcontour',0);
title('PCA forward model')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                 %
%    Source separation via GED    %
%                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generalized eigendecomposition (GED)
[evecs,evals] = eig(cov,cov);

% sort eigenvalues/vectors
[evals,sidx] = sort(diag(evals),'descend');
evecs = evecs(:,sidx);



% plot the eigenspectrum
figure(4), clf
subplot(231)
plot(evals./max(evals),'s-','markersize',15,'markerfacecolor','k')
axis square
set(gca,'xlim',[0 20.5])
title('GED eigenvalues')
xlabel('Component number'), ylabel('Power ratio (norm-\lambda)')

% component time series is eigenvector as spatial filter for data 
% (same computation as with PCA)
comp_ts = 

%% plot for comparison

% normalize time series (for visualization)
dipl_ts = dipole_data(diploc,:) / norm(dipole_data(diploc,:));
comp_ts = comp_ts / norm(comp_ts);
chan_ts = EEG.data(31,:) / norm(EEG.data(31,:));


% plot the time series
subplot(212), hold on
plot(EEG.times,.3+dipl_ts,'linew',2)
plot(EEG.times,.15+chan_ts)
plot(EEG.times,comp_ts)
legend({'Truth';'EEG channel';'GED time series'})
set(gca,'ytick',[])
xlabel('Time (a.u.)')


%% spatial filter forward model

% The filter forward model is what the source "sees" when it looks through the
% electrodes. It is obtained by passing the covariance matrix through the filter.
filt_topo = 

% Eigenvector sign uncertainty can cause a sign-flip, which is corrected for by 
% forcing the largest-magnitude projection electrode to be positive.
[~,se] = max(abs( filt_topo ));
filt_topo = filt_topo * sign(filt_topo(se));


% plot the maps
subplot(232)
topoplotIndie(lf.GainN(:,diploc), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
title('Truth topomap')

subplot(233)
topoplotIndie(filt_topo,EEG.chanlocs,'electrodes','numbers','numcontour',0);
title('GED forward model')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                 %
%    Source separation via ICA    %
%                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% run ICA and compute IC time series and maps
icvecs = jader(EEG.data,20);
ICs    = icvecs(1,:)*EEG.data;
icmaps = pinv(icvecs');
ICenergy = sum(icmaps.^2,2);

figure(5), clf

% plot component energy
subplot(231)
plot(ICenergy./max(ICenergy),'s-','markersize',15,'markerfacecolor','k')
axis square
set(gca,'xlim',[0 20.5])
title('IC energy')
xlabel('Component number'), ylabel('Power ratio (norm-\lambda)')



% plot the maps
subplot(232)
topoplotIndie(lf.GainN(:,diploc), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
title('Truth topomap')

subplot(233)
topoplotIndie(icmaps(1,:),EEG.chanlocs,'electrodes','numbers','numcontour',0);
title('IC forward model')


% plot the time series
subplot(212), hold on
plot(EEG.times,.3+dipl_ts,'linew',2)
plot(EEG.times,.15+chan_ts)
plot(EEG.times,ICs)
legend({'Truth';'EEG channel';'IC time series'})
set(gca,'ytick',[])
xlabel('Time (a.u.)')

%%


%%
% 
%   Now for real data: mouse v1
%

% a clear MATLAB workspace is a clear mental workspace, blah blah
close all; clear all

%% Source separation for stimulus on/off

% The recording is from 16 channels in mouse V1 (1=superficial, 16=deep).
% The main stimulus appeared at .5 seconds and disappeared at 1 second.

% First, load the data and visualize the ERPs on all channels

load v1_laminar.mat
csd = double(csd); % greater accuracy for double-precision data

figure(6), clf, hold on
plot(timevec,bsxfun(@plus,mean(csd,3),(1:16)'*500),'linew',2)
set(gca,'xlim',timevec([1 end]),'ytick',[])
xlabel('Time (s)'), ylabel('Channels')

% vertical lines for stim on/offsets
plot(repmat([0 .5 1],2,1),repmat(get(gca,'ylim')',1,3),'k--')




% Select two time windows that capture the stim-on and stim-off
% periods. The time windows should be short enough to be specific, but long
% enough to make robust covariance matrices.
timewin_1 = [  ]; % time in seconds
timewin_2 = [  ]; % time in seconds



% draw patches around selected time windows
h = patch(timewin_1([1 1 2 2]),[get(gca,'ylim') fliplr(get(gca,'ylim'))],'g');
set(h,'FaceAlpha',.25,'edgecolor','g')

h = patch(timewin_2([1 1 2 2]),[get(gca,'ylim') fliplr(get(gca,'ylim'))],'r');
set(h,'FaceAlpha',.25,'edgecolor','r')


%% create covariance matrices

% convert time in seconds to indices
tidxS = dsearchn(timevec',timewin_1');
tidxR = dsearchn(timevec',timewin_2');


% Compute covariance matrices over all trials (not the ERP) in those
% two time windows.
[S,R] = deal( zeros(size(csd,1)) );

for triali=1:size(csd,3)
    
    %%% S
    % extract data from all channels, the specific time range (tidxS), and this trial
    tmp = csd
    % mean-center
    tmp = bsxfun(@minus,tmp,mean(tmp,2));
    % add to covariance matrix
    S = S + tmp*tmp' / length(tmp);
    
    %%% repeat for R
    % extract data from all channels, the specific time range (tidxR), and this trial
    tmp = csd
    % mean-center
    tmp = bsxfun(@minus,tmp,mean(tmp,2));
    % add to covariance matrix
    R = R + tmp*tmp' / length(tmp);
end


% plotting
figure(7), clf
subplot(121)
imagesc(S)
axis square
set(gca,'clim',[-.5 .5]*max(abs(S(:))))
title([ 'S matrix (' num2str(timewin_1(1)) '-' num2str(timewin_1(2)) ' s.)' ])


subplot(122)
imagesc(R)
axis square
set(gca,'clim',[-.5 .5]*max(abs(S(:))))
title([ 'R matrix (' num2str(timewin_2(1)) '-' num2str(timewin_2(2)) ' s.)' ])


%% GED source separation on stim-on vs. stim-off


% GED
[V,D] = eig(S,R);

% sort eigenvals/vecs
[d,sidx] = sort(diag(D),'descend');
V = V(:,sidx);


% compute the component time series
compts = V(:,1)'*reshape(csd,size(csd,1),[]);
compts = reshape(compts,size(csd,2),size(csd,3));

% filter forward model
filt_topo = V(:,1)'*S;


%% plot the results


figure(8), clf

% eigenspectrum
subplot(221)
plot(d,'s-','markerfacecolor','k','markersize',10)
set(gca,'xlim',[0 size(csd,1)+1])
xlabel('Components'), ylabel('SNR ratio (\lambda)')


% "topography"
subplot(222), hold on
plot(filt_topo,1:16,'s-','markerfacecolor','r','markersize',10)
plot([0 0],get(gca,'ylim'),'k--')
plot(get(gca,'xlim'),[7 7],'k:')
xlabel('Component projection'), ylabel('Channel --> superficial')
set(gca,'YDir','reverse')


% time series
subplot(212), hold on
compERP = zscore(mean(compts,2));
chanERP = zscore(mean(mean(csd(6:7,:,:),3),1));

plot(timevec,compERP, timevec,chanERP,'linew',2)
set(gca,'xlim',timevec([1 end]))

% vertical lines for stim on/offsets
plot(repmat([0 .5 1],2,1),repmat(get(gca,'ylim')',1,3),'k--')
legend({'Component','Chans 6/7'})
title('Component and channel time series')
xlabel('Time (s)'), ylabel('Activity (a.u.)')

%% done.
