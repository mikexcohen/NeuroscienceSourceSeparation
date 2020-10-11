%%
%     COURSE: Advanced neuroscience techniques
%    LECTURE: Multivariate cross-frequency coupling, part 2
% Instructor: mikexcohen.com
%
%%


%% preliminaries

close all; clear

% mat file containing EEG, leadfield and channel locations
load emptyEEG
origEEG = EEG;

% indices of dipole locations
thetaDipole =  94;
gam1Dipole  = 109;
gam2Dipole  = 111;

whichOri = 1; % 1 for "EEG" and 2 for "MEG"

% plot dipole projections
figure(1), clf
clim = [-45 45];

subplot(231), topoplotIndie(-lf.Gain(:,whichOri,thetaDipole),EEG.chanlocs,'maplimits',clim,'numcontour',0,'electrodes','off','shading','interp');
title('Theta ground truth')

subplot(232), topoplotIndie(-lf.Gain(:,whichOri,gam1Dipole), EEG.chanlocs,'maplimits',clim,'numcontour',0,'electrodes','off','shading','interp');
title('Gamma/theta ground truth')

subplot(233), topoplotIndie(-lf.Gain(:,whichOri,gam2Dipole), EEG.chanlocs,'maplimits',clim,'numcontour',0,'electrodes','off','shading','interp');
title('Gamma/notheta ground truth')

%% create correlated 1/f noise time series

% correlation matrix
cormat = rand(size(lf.Gain,3));
cormat = cormat*cormat';
% cormat = cormat-min(cormat(:));
cormat = .8*( cormat./max(cormat(:)) );
cormat(1:size(lf.Gain,3)+1:end) = 1;

% eigdecomp and create correlated random data
[evecs,evals] = eig(cormat);

% 1/f random data
ps   = bsxfun(@times, exp(1i*2*pi*rand(size(lf.Gain,3),floor(EEG.pnts/2))) , .1+exp(-(1:floor(EEG.pnts/2))/200) );
ps   = [ps zeros(size(lf.Gain,3),1) ps(:,end:-1:1)];
data = 100* real(ifft(ps,[],2))' * (evecs*sqrt(evals))';

%% replace one dipole with nonstationary theta oscillation

thetafreq = 6;

% create data time series
ampl1     = 5+10*filterFGx(randn(1,EEG.pnts),EEG.srate,1,30);
freqmod1  = detrend(10*filterFGx(randn(1,EEG.pnts),EEG.srate,1,30));
thetawave = ampl1 .* sin(2*pi.*thetafreq.*EEG.times + 2*pi/EEG.srate*cumsum(freqmod1));

%% simulate gamma bursts locked (or not) to theta troughs

gammafreq = 40; % Hz

mod1 = .1+.9*(1+real(exp(1i*angle(hilbert(thetawave)))))/2;
gamma1wave = mod1 .* sin(2*pi*gammafreq*EEG.times);

ampl2      = 2+5*filterFGx(randn(1,EEG.pnts),EEG.srate,1,30);
freqmod2   = detrend(10*filterFGx(randn(1,EEG.pnts),EEG.srate,1,30));
gamma2wave = ampl2 .* sin(2*pi.*50.*EEG.times + 2*pi/EEG.srate*cumsum(freqmod2));

%% replace key dipole time series with simulated data

data(:,thetaDipole) = thetawave*10;
data(:,gam1Dipole)  = gamma1wave*10;
data(:,gam2Dipole)  = gamma2wave*20;

% project dipole data to scalp
EEG.data = detrend( data*squeeze(lf.Gain(:,whichOri,:))' )';

%% GED to identify theta component

% theta covariance
thetafilt = filterFGx(EEG.data,EEG.srate,thetafreq,5);
thetafilt = bsxfun(@minus,thetafilt,mean(thetafilt,2));
thetacov  = (thetafilt*thetafilt')/EEG.pnts;

% broadband covariance
tmpdat = bsxfun(@minus,EEG.data,mean(EEG.data,2));
bbcov  = (tmpdat*tmpdat')/EEG.pnts;

% GED
[evecsT,evals] = eig(thetacov,bbcov);

% find best component and compute filter projection
[~,maxcomp] = sort(diag(evals));
thetamap    = thetacov*evecsT;
thetamap    = thetamap(:,maxcomp(end));

% fix sign of map (max is positive)
[~,maxe] = max(abs(thetamap));
thetamap = thetamap * sign(thetamap(maxe));

% theta time series component
thetacomp = thetafilt' * evecsT(:,maxcomp(end));

% fix sign of time series according to sign of correlation with EEG
thetacomp = thetacomp * sign(corr(thetacomp,filterFGx(EEG.data(maxe,:),EEG.srate,thetafreq,9)'));

%% identify troughs and get surrounding covariance matrices

nwin = ceil(EEG.srate/thetafreq/8); % window size is 1/4 cycle (1/8 of either side)

% find troughs
troughs = find(diff(sign(diff( thetacomp )))>0)+1;
troughs(troughs<nwin+1) = [];
troughs(troughs>EEG.pnts-nwin-1) = [];


covT = zeros(EEG.nbchan);

% trough-locked covariance
for ti=1:length(troughs)
    tmpdat = EEG.data(:,troughs(ti)-nwin:troughs(ti)+nwin);
    tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
    covT   = covT + (tmpdat*tmpdat')/nwin;
end
covT = covT./ti;

%% GED to get gamma peak/trough networks

[evecs,evals] = eig(covT,bbcov);
[~,compidx]   = sort(diag(evals)); % max component

maps    = covT*evecs; % forward model of filter
gamnet1 = maps(:,compidx(end));

% fix sign
[~,idx] = max(abs(gamnet1));
gamnet1 = gamnet1 * sign(gamnet1(idx));

% plot component maps for comparison with dipole projections
figure(1)
subplot(234), topoplotIndie(thetamap, EEG.chanlocs,'numcontour',0,'electrodes','off','shading','interp');
title('Theta GED comp.')

subplot(235), topoplotIndie(gamnet1,  EEG.chanlocs,'numcontour',0,'electrodes','off','shading','interp');
title('Gamma/theta comp.')

%% get time course and reconstruct topography and sources

frex    = linspace(10,190,70);
mvarpac = zeros(size(frex));

for fi=1:length(frex)
    
    % bandpass filter trough component
    gam1comp = filterFGx(EEG.data,EEG.srate,frex(fi),15)' * evecs(:,compidx(end));
    
    % find peaks and troughs and get distances
    troughs  = find(diff(sign(diff( gam1comp )))>0)+1;
    peeks    = find(diff(sign(diff( gam1comp )))<0)+1;
    mvarpac(fi) = mean(gam1comp(peeks)) - mean(gam1comp(troughs));
end

%% more plotting

figure(2), clf

% vector of frequencies
hz = linspace(0,EEG.srate,EEG.pnts);


% plot components power spectra
subplot(331)
plot(hz,abs(fft(EEG.data'*evecsT(:,end))/EEG.pnts)), hold on
plot(hz,abs(fft(EEG.data'*evecs(:,compidx(end)))/EEG.pnts))
set(gca,'xlim',[0 90],'xtick',10:20:90), axis square
xlabel('Frequency (Hz)'), ylabel('Power (a.u.)')
legend({'Theta comp','Trough comp'})
title('Components spectra')

% plot channel power spectra (hard-coded for POz)
subplot(332)
plot(hz,abs(fft(EEG.data(30,:))/EEG.pnts),'k')
set(gca,'xlim',[0 90],'xtick',10:20:90), axis square
xlabel('Frequency (Hz)'), ylabel('Power (a.u.)')
title('Channel 30 spectrum')

% plot CFC modulation
subplot(333)
plot(frex,mvarpac,'k','linew',2)
set(gca,'xlim',[0 90],'xtick',10:20:90), axis square
xlabel('Frequency (Hz)'), ylabel('CFC strength')

% sample time courses
subplot(312)
plot(EEG.times,zscore(thetacomp), EEG.times,zscore(abs(hilbert(filterFGx(EEG.data,EEG.srate,40,25)' * evecs(:,compidx(end))))),'linew',2);
set(gca,'xlim',[2 4])
legend({'Theta component';'Gamma comp. power'})
xlabel('Time (s)'), ylabel('EEG activity (a.u.)')

subplot(313)
plot(EEG.times,zscore(filterFGx(EEG.data(30,:),EEG.srate,thetafreq,5)'), EEG.times,zscore(abs(hilbert(filterFGx(EEG.data(30,:),EEG.srate,40,25)'))),'linew',1);
set(gca,'xlim',[2 4],'ylim',[-2.5 3.5])
legend({'theta filtered';'gamma filtered'})

%% now for normal Euler-based univariate PAC (note: takes a while to run for all channels!)

pacz = zeros(EEG.nbchan,length(frex));

nperm     = 1000;
cutpoints = randsample(10:EEG.pnts-10,nperm);
permpac   = zeros(1,nperm);

% can also run only for POz (channel 30)
for chani=30%1:EEG.nbchan
    
    % filter for 6 Hz
    phase = angle(hilbert(filterFGx(EEG.data(chani,:),EEG.srate,thetafreq,5)));
    
    % filter at a higher frequency
    for fi=1:length(frex)
        pow = abs(hilbert(filterFGx(EEG.data(chani,:),EEG.srate,frex(fi),5)));
        
        % zPACd
        obspac = abs(mean( pow.*( exp(1i*phase) ) ));
        for permi=1:nperm
            permpac(permi) = abs(mean( pow([cutpoints(permi):end 1:cutpoints(permi)-1]).*( exp(1i*phase) ) ));
        end
        
        pacz(chani,fi) = (obspac-mean(permpac))/std(permpac);
    end
end


% now plot
figure(3), clf
clim = [-2.5 2.5];

subplot(221)
imagesc(1:64,frex,pacz')
xlabel('Channel'), ylabel('Frequency (Hz)')
set(gca,'clim',clim)

subplot(222), topoplotIndie(pacz(:,dsearchn(frex',40)),EEG.chanlocs,'maplimits',clim);
subplot(223), topoplotIndie(pacz(:,dsearchn(frex',65)),EEG.chanlocs,'maplimits',clim);

subplot(224)
plot(frex,pacz(30,:),'k','linew',2)
set(gca,'xlim',[0 90],'xtick',10:20:90), axis square
xlabel('Frequency (Hz)'), ylabel('CFC strength')

%% end.
