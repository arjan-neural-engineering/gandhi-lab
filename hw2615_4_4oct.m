clear all
clear space

% Arjan Singh Puniani
% BIOE 2615
% HW 4

%% Part I: Q1
% Filter the data using a high-pass Butterworth filter, use an appropriate order and
% frequency cut-off. Use the butter and filter commands. Calculate the mean and standard
% deviations of the filtered data. Justify your selection.

load 1000702VisuallyEvoked.mat

% Butterworth filter
order = 4;  cutoff = 500;
%iOrder = order/2; % only applies to bandpass
[bw_b,bw_a]=butter(order,cutoff/(0.5*samprate),'high');

channels=numel(Wideband_data);
endSample = min(cellfun(@numel,Wideband_data)); % great general function from Prof TK's HW3 sol
tVec = [1:endSample]/samprate; % time in seconds (since samprate is given in Hz)
data=zeros(endSample,channels); % (rows = samples) x (cols = ch)

for ch=1:channels
  data(:,ch)=filter(bw_b,bw_a,Wideband_data{ch}(:)); 
end

filt_avg = 1e6*mean(data,1); filt_std = 1e6*std(data,1);
% Show first [sec] seconds for channel #[tmpch]
sec = 0.25; %specify here
tmpch=7;
figure(1), clf,
% raw
subplot(311), plot(1e3*tVec(1:round(sec*samprate)),Wideband_data{tmpch}(1:round(sec*samprate))*1e3)
xlabel('Time (ms)'), ylabel('Amplitude (mV)'), title("raw data for ch #" + tmpch),
axis('tight'), grid('on')
% filtered
subplot(312), plot(1e3*tVec(1:round(sec*samprate)),data(1:round(sec*samprate),4)*1e3)
xlabel('Time (ms)'), ylabel('Amplitude (mV)'), title("high-pass filtered data for ch #" + tmpch),
axis('tight'), grid('on'),
% superimposed
subplot(313), plot(1e3*tVec(1:round(sec*samprate)),Wideband_data{tmpch}(1:round(sec*samprate))*1e3), hold on
plot(1e3*tVec(1:round(sec*samprate)),data(1:round(sec*samprate),4)*1e3), hold off
xlabel('Time (ms)'), ylabel('Amplitude (mV)'), title("4th-order high-pass filtered data atop raw for ch #" + tmpch),
axis('tight'), grid('on'), legend('raw','filtered')

figure(2)
plot(tVec,1e6*Wideband_data{tmpch}), hold on
plot(tVec,1e6*data(:,tmpch)), hold off
xlabel('Time (s)'), ylabel('Amplitude (uV)'), title("4th-order high-pass filtered data (red) superimposed atop raw data (blue) for ch #" + tmpch),
axis('tight'), grid('on'), yline(filt_avg(tmpch),'--','filtered \mu','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle'),
yline(filt_avg(tmpch)+filt_std(tmpch),'--w','+ \sigma','LabelVerticalAlignment','middle'), yline(filt_avg(tmpch)-filt_std(tmpch),'--w','- \sigma','LabelVerticalAlignment','middle'),
yline(filt_avg(tmpch)-(3*filt_std(tmpch)),'--w','- 3\sigma','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle'),
filt_avg(tmpch);
filt_std(tmpch);
%
%% Q2 & Q3
% Use the TrigON variable to find the spacing between stimuli, then re-organize the data into a new variable that has 
% dimensions of (trial_data,stimulus#). Since the trigger indicates the time of stimulus, make sure you grab 
% some time prior to stimulation (e.g. 0.5 sec below).
n_stimulations = numel(TrigON);
stim_time = mean(diff(TrigON)); stim_duration = max(TrigOFF-TrigON);
new_indeces=[1:round(samprate*stim_time)]-round(0.5*samprate); % since the samprate tells us how many samples/sec
datafilter=cell(1,channels); % Try Prof. T.K.'s way
for ch = 1:channels
   datafilter{ch} = filter(bw_b,bw_a,double(Wideband_data{ch}(:)),[],2);
end
std_mult = 3.3;
raster=zeros(size(data));
for ch=1:channels
   datafilterch = datafilter{ch};
   negThresh = mean(datafilterch) - (std_mult*(std(datafilterch)));
   indexOfCrossings = find(datafilterch < negThresh);
   id=find(diff(indexOfCrossings) > 1);
   indexOfCrossings = [indexOfCrossings(1);indexOfCrossings(id+1)];
   raster(indexOfCrossings,ch)=1;
end
psth_dt = 0.05; % should be between 25 and 200 ms
time_psth = [-0.5:psth_dt:stim_duration];
t=([0:length(new_indeces)-1]/samprate)+time_psth(1); % raster time but we want it to offset by the -0.5 delay from tVec_psth
megaRaster = zeros(numel(new_indeces),n_stimulations,channels); %newdata=zeros(numel(new_indeces),n_stimulations,channels);
for m=1:n_stimulations % find the first time t larger than the stimulus time
    stimloc(m) = find(tVec>TrigON(m),1);    % newdata(:,m)=data(stimloc(m)+new_indeces);
end
for ch = 1:channels
    for m=1:n_stimulations
        megaRaster(:,m,ch)=raster(stimloc(m)+new_indeces,ch);
    end
end

figure(3)
newtmp = 7;
plot_raster(t,megaRaster(:,1:6:n_stimulations-4,newtmp)), title("Raster for ch #" + newtmp)

figure(4)
for ch=1:channels
    subplot(sqrt(channels),sqrt(channels),ch)
    plot_raster(t,megaRaster(:,:,ch)), title("ch #" + ch),
end

%% Q4: PSTH

for m=1:length(time_psth),
    search_index=find( (t>=time_psth(m)) & (t<time_psth(m)+psth_dt) );
    psth(m,:,:)=sum(megaRaster(search_index,:,:),1);
end;

% experiment with squeeze to eliminate "1-dimensionals" in matrices
firstAvg = squeeze(mean(psth,2)); % 31 (time_psth) x 16
secondAvg = mean(firstAvg,2)/psth_dt';
figure(5)
bar(secondAvg,psth_dt), title('histogram of trial-averaged and channel-averaged  # of spikes'),
xlabel('Spikes'), ylabel('Frequency')

figure(6)
plot(time_psth,secondAvg), title('mean PSTH across all channels'),
xlabel('Time (s)'), ylabel('# spikes/bin'), xline(0,'-.','stimulus ON','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle'),

figure(7)
plot_psth(t,megaRaster(:,:,newtmp),psth_dt), title('PSTH across all channels using f(x) = plot raster'),
xlabel('Time (s)'), ylabel('# spikes/bin'), xline(0,'-.','stimulus ON','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle'),
legend('ch1','ch2','ch3','ch4','ch5','ch6','ch7','ch8','ch9','ch10','ch11','ch12','ch13','ch14','ch15','ch16','stimulus'),
