%%
%     COURSE: Advanced neuroscience techniques
%    LECTURE: Multivariate cross-frequency coupling
% Instructor: mikexcohen.com
%
%%


%%
%
% Example PAC in human nucleus accumbens recording
%
%


% eight seconds of data recorded from 
% the human nucleus accumbens, sampled at 1 kHz
load accumbens_eeg.mat
srate = 1000;
npnts = length(eeg);

% first let's see the data
figure(1)
plot(eeg)

%% plot PAC

% specify frequency bands for phase and amplitude
freq4phase = 10; % in Hz
freq4power = 70; 

% get phase values
phasefilt = filterFGx(eeg,srate,freq4phase,.5);
phase = angle(hilbert(phasefilt));

% get power values (note: 'power' is a built-in function so we'll name this variable 'amp')
ampfilt = filterFGx(eeg,srate,freq4power,40);
amp = abs(hilbert(ampfilt)).^2;

% plot power and phase
figure(2), clf
plot(phase)
hold on
plot((amp-mean(amp))/std(amp),'r')
% hint: zoom in on part of these data!

% plot power as a function of phase
figure(3), clf
polarplot(phase,amp,'.')

%% alternative visualization as histogram

n_hist_bins = 30;

phase_edges = linspace(min(phase),max(phase),n_hist_bins+1);
amp_by_phases = zeros(1,n_hist_bins);

for i=1:n_hist_bins-1
    amp_by_phases(i) = mean(amp(phase>phase_edges(i) & phase<phase_edges(i+1)));
end

figure(4)
bar(phase_edges(1:end-1),amp_by_phases,'histc');
set(gca,'xlim',[phase_edges(1) phase_edges(end)])
xlabel([ 'Phase at ' num2str(freq4phase) ' Hz (rad.)' ])
ylabel([ 'Power at ' num2str(freq4power) ' Hz' ])


%% significance through permutation testing 

% compute amplitude-modulated phase
observed_modulation = abs(mean(amp.*exp(1i*phase)));

% now run permutation testing to see the significance of this interaction
num_iterations = 2000;
permuted_results = zeros(1,num_iterations);


for ni=1:num_iterations
    % randomly resort power values
    cutpoint = randsample(round(npnts/10):round(npnts*.9),1);
    permuted_results(ni) = abs(mean(amp([ cutpoint:end 1:cutpoint-1 ]).*exp(1i*phase)));
end

% print value
z_modulation = (observed_modulation-mean(permuted_results))/std(permuted_results);
disp([ 'Z-modulation of ' num2str(freq4power) ' Hz power modulated by ' num2str(freq4phase) ' Hz phase is ' num2str(z_modulation) '.' ])

% show histogram of permuted and real values
figure(5), clf
histogram(permuted_results,50);
hold on
plot([observed_modulation observed_modulation],get(gca,'ylim')/2,'m-p','linewi',4,'markersize',16);
legend({'histogram of permuted values';'observed value'})
xlabel('Absolute modulate (arb. units depending on power values)')
ylabel('Count')


%% generate a full 2D plot of CFC

% define frequencies for phase and for amplitude
phas_freqs =  5:2:20;
ampl_freqs = 30:5:150;

% number of iterations used for permutation testing
n_iter = 200; 

% initialize output phase-amplitude matrix
phaseamp = zeros(length(phas_freqs),length(ampl_freqs));


% loop over frequencies for phase
for lower_fi=1:length(phas_freqs)
    
    % get phase values
    phasefilt = filterFGx(eeg,srate,phas_freqs(lower_fi),phas_freqs(lower_fi)*.4);
    phase = angle(hilbert(phasefilt));
    
    for upper_fi=1:length(ampl_freqs)
        
        % get power values (note: 'power' is a built-in function so we'll name this variable 'amp')
        ampfilt = filterFGx(eeg,srate,ampl_freqs(upper_fi),ampl_freqs(upper_fi)*.78);
        amplit = abs(hilbert(ampfilt)).^2;
        
        % calculate observed modulation index
        modidx = abs(mean(amplit.*exp(1i*phase)));

        % now use permutation testing to get Z-value
        bm = zeros(1,length(n_iter));
        for bi=1:n_iter
            cutpoint = randsample(round(npnts/10):round(npnts*.9),1);
            bm(bi) = abs(mean(amplit([ cutpoint:end 1:cutpoint-1 ]).*exp(1i*phase)));
        end

        % the value we use is the normalized distance away from the mean of
        % boot-strapped values
        phaseamp(lower_fi,upper_fi) = (modidx-mean(bm))/std(bm);
    end % end upper frequency loop (for amplitude)
end % end lower frequency loop (for phase)


% plot it! let's try a contour map
figure(6), clf
contourf(phas_freqs,ampl_freqs,phaseamp',40,'linecolor','none')
set(gca,'clim',[-3 3])
xlabel('Frequency for phase')
ylabel('frequency for amplitude')
title('Map of phase-amplitude coupling in human nucleus accumbens')
% note: red colors mean more coupling

%% 

close all; clear

%%
%
% Positive PAC from non-sinusoidal waveforms
%
%


%% signal parameters

srate = 1000;
N = srate*2;
time = (0:N-1)/srate;

lowerfreq = 6;
upperfreq = 40; % Hz

%% create the signal

% option 1: (not CFC) very simple
signal = sin(2*pi*lowerfreq*time) + randn(1,N)/100;
signal = signal + (signal>.5) .* sin(2*pi*upperfreq*time)/3;


% option 2: (with CFC) more physiological
ampl1     = 5+10*filterFGx(randn(1,N),srate,1,30);
freqmod1  = detrend(10*filterFGx(randn(1,N),srate,1,30));
thetawave = ampl1 .* sin(2*pi*lowerfreq*time + 2*pi/srate*cumsum(freqmod1));
mod1      = .1+.9*(1+real(exp(1i*angle(hilbert(thetawave)))))/2;
gammawave = mod1 .* sin(2*pi*upperfreq*time);
signal    = thetawave + gammawave*2;


% 
% option 3: (not CFC) bit of a square at the top
signal = sin(2*pi*lowerfreq*time);
signal(signal>.8) = 1;


%%

%%% Compute PAC
phase_frex = linspace(2,15,20);
power_frex = linspace(20,80,60);

pac = zeros(length(power_frex),length(phase_frex));

for lowfi=1:length(phase_frex)
    
    % get phase time series
    phasets = exp(1i*angle(hilbert(filterFGx(signal,srate,phase_frex(lowfi),1))));
    
    % loop over amplitude frequencies
    for hifi=1:length(power_frex)
        
        % get power time series
        powts = abs(hilbert(filterFGx(signal,srate,power_frex(hifi),lowerfreq*2))).^2;
        
        % PAC
        pac(hifi,lowfi) = abs(mean( powts.*phasets ));
    end
end


%%% visualize

figure(1), clf
subplot(211)
plot(time,signal,'linew',3)
xlabel('Time (s)')
box off, set(gca,'ytick',[])

subplot(212)
contourf(phase_frex,power_frex,pac,40,'linecolor','none')
xlabel('Phase freq. (Hz)')
ylabel('Ampl. freq. (Hz)')


%% 

close all; clear

%%
%
% Cross-frequency coupling via GED
%
%

% go to other script...

%% done.
