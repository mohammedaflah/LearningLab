%%
%     COURSE: Solved problems in neural time series analysis
%    SECTION: Time-domain analyses
%      VIDEO: Project 2-2: ERP peak latency topoplot
% Instructor: sincxpress.com
%
%%

%% Loop through each channel and find the peak time of the ERP between 100 and 400 ms. 
%   Store these peak times in a separate variable, and then make a
%   topographical plot of the peak times. 
load sampleEEGdata.mat

% Define time boundaries and convert to indices
timeboundaries = [100 400];
timeidx = dsearchn(EEG.times',timeboundaries');

[~,maxERPtime] = max( mean(EEG.data (:,timeidx(1):timeidx(2),:), 2));

%Convert indices back to indices
maxERPtime = EEG.times(maxERPtime+timeidx(1)-1);

% make plot
figure(1), clf
subplot(121)
topoplotIndie(maxERPtime,EEG.chanlocs,'numcontour',4,'electordes','numbers');
title({'Unfiltered ERP peak time';[' (' num2str(timeboundaries(1)) '-' num2str(timeboundaries(2)) 'ms)']}); 
set(gca,'clim',timeboundaries)
colormap hot
colorbar

%% Repeat for a low-pass filtered ERP.


lowcut = 15;
filttime = -.3:1/EEG.srate:.3;
filtkern = sin(2*pi*lowcut*filttime) ./ filttime;

% adjust NaN and normalize filter to unit-gain
filtkern(~isfinite(filtkern)) = max(filtkern);
filtkern = filtkern./sum(filtkern);

% windowed sinc filter
filtkern = filtkern .* hann(length(filttime))';

%
erp = zeros(EEG.nbchan,EEG.pnts);
for chani=1:EEG.nbchan
    erp(chani,:) = filtfilt(filtkern,1,double(mean(EEG.data(chani,:,:),3)));
end

[~,maxERPtime] = max( abs(erp(:,timeidx(1):timeidx(2),:), 2));

%indices to ms
maxERPtime = EEG.times(maxERPtime+timeidx(1)-1);

% make plot
subplot(122)
topoplotIndie(maxERPtime,EEG.chanlocs,'numcontour',4,'electordes','numbers');
title({'filtered ERP peak time';[' (' num2str(timeboundaries(1)) '-' num2str(timeboundaries(2)) 'ms)']}); 
set(gca,'clim',timeboundaries)
colormap hot
colorbar
