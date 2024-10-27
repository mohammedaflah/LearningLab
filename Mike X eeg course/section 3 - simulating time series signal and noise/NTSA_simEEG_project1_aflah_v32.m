%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Simulating EEG data
%      VIDEO: Project 1: Channel-level EEG data
% Instructor: sincxpress.com
%
%%

%%% INSTRUCTIONS:
% The goal of this assignment is to simulate time series data
% that can be used to test time-series analysis methods.
% For each section below: 
%   1) Complete the MATLAB code
%   2) Put the data into the EEG structure
%      - Make sure all relevant fields are accurate (EEG.data, EEG.pnts, EEG.trials, EEG.srate, EEG.nbchan, EEG.times)
%   3) Use function plot_simEEG to plot some data

%% 1) pure phase-locked sine wave

% parameters
EEG.srate  = 500; % sampling rate in Hz
EEG.pnts   = 1500;
EEG.trials = 30;
EEG.nbchan = 23;

sinefreq = 6.75; % in Hz

% time vector
EEG.times = (0:EEG.pnts-1)/EEG.srate;


% loop over channels and create data
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        % data as a sine wave
        EEG.data(chani,:,triali) = sin(2*pi*sinefreq*EEG.times);
    end
end


% plot an ERP from one channel
figure(1), clf
plot(EEG.times,squeeze(mean(EEG.data(10,:,:),3)),'linew',2)
xlabel('Time (s)'), ylabel('Activity')
title('ERP from channel 10')


% the function below takes at least one argument (EEG),
% and optionally a second argument (channel number),
% and optionally a third argument (figure number)
plot_simEEG(EEG,2,3)

%% 2) Non-phase-locked sine wave

% hint: copy/paste the code above but add something inside the sine 
%       function on each trial.
EEG.srate  = 500; % sampling rate in Hz
EEG.pnts   = 1500;
EEG.trials = 30;
EEG.nbchan = 23;

sinefreq = 6.75; % in Hz

% time vector
EEG.times = (0:EEG.pnts-1)/EEG.srate;


% loop over channels and create data
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        % data as a sine wave
        EEG.data(chani,:,triali) = sin(2*pi*sinefreq*EEG.times+ 2*pi*randn);
    end
end


% plot an ERP from one channel
figure(1), clf
plot(EEG.times,squeeze(mean(EEG.data(10,:,:),3)),'linew',2)
xlabel('Time (s)'), ylabel('Activity')
title('ERP from channel 10')
plot_simEEG(EEG,2,4)

%% 3) multisine waves

% list of frequencies and corresponding amplitudes
frex = [ 3 5 16 ];
amps = [ 2 4 5  ];


% loop over channels and trials
for chani=1:EEG.nbchan
    for triali=1:EEG.trials

        sinewave = zeros(1,EEG.pnts);
        for si=1:numel(frex)
            sinewave = sinewave + amps(si)* sin(2*pi*frex(si)*EEG.times);
        end
        
        % data as a sine wave plus noise
        EEG.data(chani,:,triali) = sinewave + randn ;
    end
end

plot_simEEG(EEG,2,5)

%% 4) nonstationary sine waves

% loop over channels and trials
for chani=1:EEG.nbchan
    for triali=1:EEG.trials

        % hint: instantaneous frequency via interpolated random numbers
        freqmod = 20*interp1(rand(1,10),linspace(1,10,EEG.pnts));
        signal  = sin( 2*pi * ((EEG.times + cumsum(freqmod))/EEG.srate) );

        EEG.data(chani,:,triali) = signal ;
    end
end

plot_simEEG(EEG,2,6)


%% 5) transient oscillations w/ Gaussian


peaktime = 1; % seconds
width = .12;
sinefreq = 7;

gaus = exp( -(EEG.times-peaktime).^2 / (2*width^2) );

% loop over channels and trials
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        cosw = cos(2*pi*sinefreq*EEG.times)
        % then multiply the gaussian by a sine wave
        EEG.data(chani,:,triali) = cosw .* gaus ;
    end
end

plot_simEEG(EEG,2,7)




%% 6) repeat #3 with white noise

% list of frequencies and corresponding amplitudes
frex = [ 3 5 16 ];
amps = [ 2 4 5  ];


% loop over channels and trials
for chani=1:EEG.nbchan
    for triali=1:EEG.trials

        sinewave = zeros(1,EEG.pnts);
        for si=1:numel(frex)
            sinewave = sinewave + amps(si)* sin(2*pi*frex(si)*EEG.times);
        end
        
        % data as a sine wave plus noise
        EEG.data(chani,:,triali) = sinewave + 5 * randn(size(sinewave)) ;
    end
end

plot_simEEG(EEG,2,8)




%% 7) repeat #5 with 1/f noise

noiseamp = .3;
peaktime = 1; % seconds
width = .12;
sinefreq = 7;

gaus = exp( -(EEG.times-peaktime).^2 / (2*width^2) );

% loop over channels and trials
for chani=1:EEG.nbchan
    for triali=1:EEG.trials

        cosw = cos(2*pi*sinefreq*EEG.times * 2*pi*rand);
        % then multiply the gaussian by a sine wave
        
        ed = 50;
        as = rand(1,floor(EEG.pnts/2)-1) .* exp(-(1:floor(EEG.pnts/2)-1)/ed);
        as = [as(1) as 0 as(:,end:-1:1)];
        EEG.data(chani,:,triali) = cosw .* gaus ;
        %Fourier coefficients
        fc = as .* exp(1i*2*pi*rand(size(as)));

        % Inverse Fourier transform to create the noise
        noise = real(ifft(fc)) * EEG.pnts ;

         EEG.data(chani,:,triali) = gaus .* cosw + noiseamp * noise ;
    end
end

plot_simEEG(EEG,2,9)








%%

