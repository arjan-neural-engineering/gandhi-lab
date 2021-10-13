clear all
clear space
close all

load bl_sc_072315_1_mcell_spikelfpcSC_fulltrial.mat

possibles = numel(data); % possible trials
channels = numel(data(1).spikeTimestamps);

%% user-defined parameters
% restricting time to the delay period
postTargBuffer = 200; % how many ms after target appears to begin analysis
% movmean smoothing params
smoothamt = 10; %k windows
smoothMult = 1.5; % * smoothamt
% thresholds
velocThresh = 15; % degrees per second
% findpeak parameters
minPeakH = 7.5; % (deg/sec), lower limit for peak microsaccade speed
maxPeakH = 90; % (deg/sec), upper limit for peak microsaccade speed
minPeakSep = 50; % number of time bins separating microsaccades
minProminen = 20; % PROMINENCE aka the relative importance
% turn plotting on (1) / off (0)
plotting = 0;

%% main loop through possible trials
for trial=1:possibles
    if (data(trial).inTarg)
        pks = 0; locs = 0;
        getTarg = data(trial).targtimes;
        getDelay = getTarg+data(trial).delays;
        idxrange = [getTarg+postTargBuffer:getDelay];
        getXYs = data(trial).gazePosition(:,idxrange); % at valid trial, get the x and y positions -- stacked x atop, y on bottom
        [thetas,rhos]=cart2pol(getXYs(1,:),getXYs(2,:)); % convert to polar coordinates
        velocVec = [NaN diff(rhos)*1e3]; % NaN to match dims after diff + convert ms to sec
        smoothV = abs(movmean(velocVec,smoothamt,2));
        smootherV = abs(movmean(smoothV,smoothMult*smoothamt,2)); % if smoothMult = 2, doubles the windows
        tVec = [1:numel(smootherV)];
        %% find peaks function
        [pks,locs]= findpeaks(smootherV,'MinPeakHeight',minPeakH,'MinPeakDistance',minPeakSep);
        if ~isempty(pks)        %% only plot if findpeaks doesn't return anything
            getVelocInd = find(smootherV>velocThresh);
            keepFirstIdx(trial)=locs(1); % KEEP only the first detected microsaccade
            if plotting
                figure
                findpeaks(smootherV,'MinPeakHeight',minPeakH,'MinPeakDistance',minPeakSep,'Annotate','extents'),
                plot(tVec,smootherV,'b'), axis('tight'), xlabel('time (ms)'), ylabel('pupil velocity (deg/sec)'),
                title("delay time-restricted window in trial #" + trial),
                if ~isempty(locs)
                    for i=1:numel(locs)
                        xline(locs(i),'--b','microsacc');
                    end
                end
                legend('extra smoothed (increased window size)','peak microsacc')
            end
        end
    end
end

%% correct the index so that it matches absolute time
for i=1:numel(keepFirstIdx)
    if (keepFirstIdx(i))
       disp("trial #" + i + " had a microsaccade at index " + keepFirstIdx(i)) ;
    end
end

%% JUNK
% in case you catch exceptions
try
    [pks,locs] = findpeaks(smoothV,'MinPeakHeight',minPeakH,'Annotate','extents');
catch Invalid MinPeakHeight
    try
        [pks,locs] = findpeaks(smoothV,'MinPeakProminence',minProminen,'Annotate','extents');
    catch
        disp("no microsaccade on trial #" + trial);
    end
end