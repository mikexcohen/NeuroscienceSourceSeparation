%%
%     COURSE: Advanced neuroscience techniques
%    LECTURE: Spectral source separation
% Instructor: mikexcohen.com
%
%%


%%
% Intro to "static" spectral analysis via FFT
%

%% Generate a multispectral noisy signal

% simulation parameters
srate = 1234; % in Hz
npnts = srate*2; % 2 seconds
time  = (0:npnts-1)/srate;

% frequencies to include
frex  = [ 12 18 ];


% loop over frequencies to create signal
signal = zeros(size(time));
for fi=1:length(frex)
    signal = signal + fi*sin(2*pi*frex*time);
end


% add noise to the signal
noiselevel = 0;
signal = signal  noiselevel * randn(size(signal));


% plot the time-domain signal
figure(1), clf
subplot(211)

plot(time,signal,'k')
xlabel('Time (s)'), ylabel('Amplitude')
title('Time domain')


%% spectral analysis

% amplitude spectrum via Fourier transform
signalX = fft();
signalAmp = 2*abs(signalX)/npnts;

% vector of frequencies in Hz
hz = linspace(0,srate/2,floor(npnts/2)+1);


% and plot
subplot(212)
stem(hz,signalAmp(1:length(hz)),'ks','linewidth',2,'markersize',10)
set(gca,'xlim',[0 max(frex)*3])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Frequency domain')


%%% Question: Is it easier to understand the features of the signal in the
%             time-domain or in the frequency-domain?
%
%%% Question: Which domain is more robust to noise? 
%             Try changing the "noiselevel" variable.

%%

%% example static power spectrum in real EEG data

% load data
load EEGrestingState.mat
N = length(eegdata);

% time vector
timevec = (0:N-1)/srate;

% plot the data in the time domain
figure(2), clf
plot(timevec,eegdata,'k')
xlabel('Time (seconds)'), ylabel('Voltage (\muV)')
title('Time domain')

%%% Question: Can you identify spectral features in the time domain?
%             Zoom in and count!
zoom on

%% "static" FFT over entire period

eegpow = abs( fft(eegdata)/N ).^2;
hz = linspace(0,srate/2,floor(N/2)+1);

figure(3), clf
plot(hz,eegpow(1:length(hz)),'k')
set(gca,'xlim',[0 60],'yscale','log')
xlabel('Frequency (Hz)'), ylabel('Power')

%%% Question: MXC claims there are 3 prominent features in this power spectrum.
%             What are they and do you see them in this spectrum?!

%% MATLAB pwelch function (window+FFT)

% create Hann window
winsize = 2*srate; % 2-second window
hannw = .5 - cos(2*pi*linspace(0,1,winsize))./2;

% number of FFT points (frequency resolution)
nfft = srate*100;

figure(4), clf
pwelch(eegdata,hannw,round(winsize/4),nfft,srate); % for MATLAB
%pwelch(eegdata,hannw,.5,nfft,srate); % Octave uses this line instead
set(gca,'xlim',[0 60])
title('Power via Welch''s method')

%%% Question: How does this figure compare to the previous one?

%%


%% Effects of non-stationarities on static power spectrum

% simulation details
srate = 1000;
t = 0:1/srate:10;
n = length(t);


% create signals (sine wave and linear chirp)
f  = [2 10];
ff = linspace(f(1),mean(f),n);
signal1 = (2*pi.*ff.*t); % create a sine wave
signal2 = sin(2*pi.*mean(ff).*t);

% Fourier spectra
signal1X = fft(signal1)/n;
signal2X = fft(signal2)/n;
hz = linspace(0,srate/2,floor(n/2));



% plot the signals in the time domain
figure(5), clf
subplot(211)
plot(t,signal1,'b'), hold on
plot(t,signal2,'r')
xlabel('Time (sec.)'), ylabel('Amplitude')
set(gca,'ylim',[-1.1 1.1])
title('Time domain')



% and their amplitude spectra
subplot(212)
stem(hz,2*abs(signal1X(1:length(hz))),'.-','linewidth',2), hold on
stem(hz,2*abs(signal2X(1:length(hz))),'r.-','linewidth',2)
xlabel('Frequency (Hz)'), ylabel('Amplitude')
set(gca,'xlim',[0 20])
title('Frequency domain')

%% 

%%
% Intro to time-frequency analysis
%


%% create signal (chirp) used in the following examples

% simulation details and create chirp
fs     = 1000; % sampling rate
time   = 0:1/fs:5;
npnts  = length(time);
f      = [10 30]; % frequencies in Hz
ff     = linspace(f(1),mean(f),npnts);
signal = sin(2*pi.*ff.*time);



figure(6), clf

% plot signal
subplot(511)
plot(time,signal,'k','linewidth',2)
xlabel('Time (s)'), ylabel('Amplitude')
title('Time domain signal')

% compute power spectrum
sigpow = 2*abs(fft()/npnts).^2;
hz     = linspace(0,fs/2,floor(npnts/2)+1);

% and plot
subplot(512)
plot(hz,sigpow(1:length(hz)),'k','linewidth',2)
xlabel('Frequency (Hz)'), ylabel('Power')
set(gca,'xlim',[0 80])


%% Creating a time-frequency plot using filter-Hilbert

% frequencies to extract
minfreq =  5; % Hz
maxfreq = 40; % Hz
numfrex = 50;

% Complex Morlet wavelet 
% convolution

% the vector of frequencies
frex = linspace(minfreq,maxfreq,numfrex);

% Gaussian filter width
fwhm = 6; % in Hz


% initialize output matrix
tf = zeros(numfrex,npnts);

for fi=1:numfrex
    
    % filter data
    fdat = filterFGx % use the function help file to figure out how to use this function
    
    % get power (amplitude squared)
    ampts = abs(hilbert( fdat )).^2;
    
    % put power time series into tf matrix for this frequency
    tf = ampts;
end

%
% plot!
subplot(5,1,3:5)
contourf(time,frex,tf,40,'linecolor','none')
% set(gca,'xlim',[-.1 1.2],'clim',[-1 1]*20) % optional figure settings
xlabel('Time (s)')
ylabel('Frequency (Hz)')


%% example in real V1 LFP data

load v1.mat

% plot the data in the time domain
figure(7), clf
subplot(411), hold on
plot(timevec,mean(data,2),'k','linewidth',2)
plot([0 0],get(gca,'ylim'),'k--')
plot([.5 .5],get(gca,'ylim'),'k--')
title('Time domain')
xlabel('Time (s)'), ylabel('Voltage (\muV)')
set(gca,'xlim',[-.2 1.4])


% static power spectrum
static = mean( abs(fft(data)/size(data,1)).^2 ,2);
hz = linspace(0,srate,size(data,1));
subplot(412)
plot(hz,static,'k','linewidth',2)
title('Frequency domain')
xlabel('Frequency (Hz)')
set(gca,'xlim',[0 80],'yscale','log','ylim',[1e2 1e4])

%% Creating a time-frequency plot using filter-Hilbert

clear % a clear MATLAB workspace is a clean mental workspace

% load data
load v1_laminar.mat
whos % check out the data

% only from one channel
chan2use = 7; % corresponds to layer-IV


% frequencies to extract (what are reasonable parameters to pick here?)
minfreq = 
maxfreq = 
numfrex = 

% the vector of frequencies
frex = linspace(minfreq,maxfreq,numfrex);

% Gaussian filter width
fwhm = 6; % in Hz

% baseline time window
% QUESTION: is this really a good choice of *baseline* time window?
base_timewin = [ 1 2 ];  % in seconds

bidx = dsearchn(timevec',base_timewin');


% initialize output matrix
tf = zeros(numfrex,length(timevec));

for fi=1:numfrex
    
    % filter data
    fdat = filterFGx( squeeze(csd(chan2use,:,:))' ,srate,frex(fi),fwhm)';
    
    % get power (amplitude squared)
    ampts = abs(hilbert( fdat )).^2;
    
    % the previous line gives power over time for each trial; now you need to average over trials
    ampts = 
    
    % get baseline
    basepow = mean(ampts(bidx(1):bidx(2)));
    
    % dB normalization
    tf(fi,:) = 10*log10( ampts/basepow );
end

%
% plot!
figure(11), clf
contourf(timevec,frex,tf,40,'linecolor','none')
set(gca,'xlim',[-.1 1.2],'clim',[-1 1]*20)
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% done.
