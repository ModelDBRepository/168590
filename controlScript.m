%%  Control Script

clear, close all 

n=2; %Number of cells used.      
numberTrials=4; %number of trials to average over
tfinal=10000; %length of each trial in ms 
stimulation=1;  %set to 0 if there is no sinusoidal component
frequency=4; %frequency of sinusoidal component in Hz                
T=1000/frequency; %period in ms
noPeriodsSampled=tfinal/T;   %number of periods over which PSTHs are averaged          
dt=0.05; %sampling period in ms
t=0:dt:tfinal; %time vector;  
CCGbinsize=0.5; %ms
CCGmaxlags=ceil(T/2); %absolute value of max and min lag when producing CCGs 
Rmaxlags=ceil(T/2);   %absolute value of max and min lag when calculating R coefficients 
noPSTHbins=12;   
PSTHbinsize=T/noPSTHbins; %ms  
c=0.25; 
e=1; 
trainingSession=0;   %set to 0 to avoid retraining weights.
trainingTime=1000000; %length of weight training. 


if trainingSession 
    
fixedNoiseTrain=normrnd(0,1,1, trainingTime/dt+1);   %shared noise between SP cells
randomNoiseTrain=normrnd(0,1,1,trainingTime/dt+1);    %noise that is unique for each SP cell

%NOTE: since we are only training one SP cell, the distinction between shared and unique noise is irrelevant.
%However, the above variables need to be set anyway. 

trialIndex=0;  %To distinguish this trial as a training trial, it has a trialIndex of zero. 

[~, ~, weightMatrix]=model(1, 1, frequency,stimulation, dt, trainingTime, randomNoiseTrain,fixedNoiseTrain, c, e, trainingSession, trialIndex); 

cd('/Users/bensimmonds/Documents/MATLAB/Noise Correlation Reduction Model/Weights'); 

save(['globalWeightMatrixe' num2str(e) num2str(frequency) 'Hz.mat'], 'weightMatrix'); 


end 

%End of training session (if applicable) 

%Now, beginning of data-collecting trials: 

trainingSession=0;

 
%Preallocate matrix space for PSTH matrices.   

localBinnedPSTH=zeros(noPeriodsSampled,noPSTHbins+1,n, numberTrials);
globalBinnedPSTH=zeros(noPeriodsSampled,noPSTHbins+1,n, numberTrials);

localBinnedAverageSamples=zeros(noPSTHbins+1, numberTrials);
globalBinnedAverageSamples=zeros(noPSTHbins+1, numberTrials); 


%Raster matrices

globalRaster=cell(1);
localRaster=cell(1); 
localChirpRaster=cell(1);
globalChirpRaster=cell(1); 


%Preallocate matrices for CCGs.

localCCGBinned=zeros(numberTrials, floor(tfinal/CCGbinsize)+1, n);
globalCCGBinned=zeros(numberTrials, floor(tfinal/CCGbinsize)+1, n);

%Perform following set of functions for each trial: 

for trialIndex=1:numberTrials;
   
    %Create noise train shared by both SP cells: 

    fixedNoiseTrain=normrnd(0,1,1, length(t)); 

    for cellNumber=1:n; 
    
        %Create noise train unique to cell n: 
    
        randomNoiseTrain=normrnd(0,1,1,length(t));  

        %Obtain spike trains by running model: 

        [localV, localSpikeTrain]=model(cellNumber, 0,frequency, stimulation,  dt, tfinal,  randomNoiseTrain,fixedNoiseTrain, c, e, trainingSession, trialIndex); 
        [globalV, globalSpikeTrain]=model(cellNumber, 1,frequency, stimulation,  dt, tfinal,  randomNoiseTrain,fixedNoiseTrain, c, e, trainingSession,  trialIndex); 
 
        %Bin local and global spike trains for each period of spike train to be sampled (for use with PSTHs). 
 
        for sampleNumber=1:noPeriodsSampled;
   
            highPeriodLimit=tfinal-(sampleNumber-1)*T;
            lowPeriodLimit=tfinal-T-(sampleNumber-1)*T;
    
            localLastPeriodIndices=sort(find(localSpikeTrain >lowPeriodLimit & localSpikeTrain <= highPeriodLimit)) ;
            localSpikeTrainPeriod=localSpikeTrain(localLastPeriodIndices) ;
            localBinnedPSTH(sampleNumber,:,cellNumber, trialIndex)=hist(localSpikeTrainPeriod, lowPeriodLimit:T/noPSTHbins:highPeriodLimit); 
    
            globalLastPeriodIndices=sort(find(globalSpikeTrain>lowPeriodLimit & globalSpikeTrain<=highPeriodLimit)); 
            globalSpikeTrainPeriod=globalSpikeTrain(globalLastPeriodIndices);  
            globalBinnedPSTH(sampleNumber,:,cellNumber, trialIndex)=hist(globalSpikeTrainPeriod,lowPeriodLimit:T/noPSTHbins:highPeriodLimit); 
   
        end 
  
        % Calculate average spike train over the sample periods for cell cellNumber (for use with PSTHs). 
 
        localBinnedAverageSamples(:, cellNumber,trialIndex)=sum(localBinnedPSTH(:,:,cellNumber, trialIndex),1)/noPeriodsSampled;
        globalBinnedAverageSamples(:, cellNumber, trialIndex)=sum(globalBinnedPSTH(:,:,cellNumber, trialIndex),1)/noPeriodsSampled; 
 

        %CCG Binning for this trial: 

        localCCGBinned(trialIndex,:,cellNumber)=hist(localSpikeTrain, 0:CCGbinsize:tfinal);
        globalCCGBinned(trialIndex,:,cellNumber)=hist(globalSpikeTrain, 0:CCGbinsize: tfinal);

    end 

end


%Now, generation of CCGs, PSTHs and R coefficients. 

%LOCAL 

%Calculate local CCGs for display:
 
localCCGRaw=CCGraw(localCCGBinned, CCGbinsize, CCGmaxlags);
localCCGSignal=CCGsignal(localCCGBinned, CCGbinsize, CCGmaxlags);
localCCGNoise=localCCGRaw-localCCGSignal; 

%Plot local raw versus local noise CCGs: 
 
figure;  plot(-CCGmaxlags/2:0.5:CCGmaxlags/2, localCCGRaw, -CCGmaxlags/2:0.5:CCGmaxlags/2,localCCGNoise);
title(['Local noise CCG versus local raw CCG with frequency ' int2str(frequency) 'Hz (Number of trials=' int2str(numberTrials) ')']); 
axis([-60 60 -20 60]); 
xlabel('Lag(ms)');
ylabel('Coincidences per second');
legend('Raw', 'Noise'); 
 

%Calculate local R coefficients: 

%First need to calculate local CCGs specifically for R calculation (because may involve
%different input parameters than CCGs created for display - e.g. maxlags) 

[localCCGRawR]=CCGraw(localCCGBinned, CCGbinsize, Rmaxlags); 
[localCCGSignalR]=CCGsignal(localCCGBinned, CCGbinsize, Rmaxlags); 
localCCGNoiseR=localCCGRawR-localCCGSignalR;

localACG1BinnedR=repmat(localCCGBinned(:,:,1),[1,1,2]);
localACG2BinnedR=repmat(localCCGBinned(:,:,2),[1,1,2]);

[localACG1RawR]=CCGraw(localACG1BinnedR, CCGbinsize, Rmaxlags);
[localACG2RawR]=CCGraw(localACG2BinnedR,  CCGbinsize, Rmaxlags);

[localACG1SignalR]=CCGsignal(localACG1BinnedR,  CCGbinsize, Rmaxlags);
[localACG2SignalR]=CCGsignal(localACG2BinnedR, CCGbinsize, Rmaxlags);

localACG1NoiseR=localACG1RawR-localACG1SignalR;
localACG2NoiseR=localACG2RawR-localACG2SignalR;

localRRaw=crossCorrelationCoefficient(localCCGRawR, localACG1RawR, localACG2RawR)

localRSignal=crossCorrelationCoefficient(localCCGSignalR, localACG1SignalR, localACG2SignalR) 

localRNoise=crossCorrelationCoefficient(localCCGNoiseR, localACG1NoiseR, localACG2NoiseR)


%Plot local PSTH (are plotting PSTH only for cell 1):   
 
figure; plot(0:PSTHbinsize:T, mean(localBinnedAverageSamples(:,1,:),3)/(PSTHbinsize/1000), '*') ;

title('Local Input PSTH for Sample Cell');
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
axis([0 T 0 100]) ;


% Plot sample local voltage trace from last trial (ie. from trial=trialNumber)

figure; plot(t, localV);
xlabel('Time(ms)');
ylabel('Voltage(mV)');
axis([0 4*T -80 -63]);
title(['Voltage trace for superficial cell under local input from a sample trial (trial number ' int2str(numberTrials) ')']); 


%GLOBAL

%Calculate global CCGs for display:

[globalCCGRaw]=CCGraw(globalCCGBinned,  CCGbinsize, CCGmaxlags);
[globalCCGSignal]=CCGsignal(globalCCGBinned, CCGbinsize, CCGmaxlags);
globalCCGNoise=globalCCGRaw-globalCCGSignal; 


%Plot global raw versus global noise CCGs: 
 
 figure;  plot(-CCGmaxlags/2:0.5:CCGmaxlags/2, globalCCGRaw, -CCGmaxlags/2:0.5:CCGmaxlags/2,globalCCGNoise);
 title(['Global noise CCG versus global raw CCG with frequency ' int2str(frequency) 'Hz (Number of trials=' int2str(numberTrials) ')' ]); 
 xlabel('lag(ms)');
  axis([-60 60 -10 30]); 
 ylabel('Coincidences per second');
 legend('Raw', 'noise'); 


% Global R calculations:

%First need to calculate global CCGs specifically for R calculation (because may involve
%different input parameters than CCGs created for display-e.g. maxlags) 

globalCCGRawR=CCGraw(globalCCGBinned, CCGbinsize, Rmaxlags); 
globalCCGSignalR=CCGsignal(globalCCGBinned, CCGbinsize, Rmaxlags); 
globalCCGNoiseR=globalCCGRawR-globalCCGSignalR;

globalACG1BinnedR=repmat(globalCCGBinned(:,:,1),[1,1,2]);
globalACG2BinnedR=repmat(globalCCGBinned(:,:,2),[1,1,2]);

globalACG1RawR=CCGraw(globalACG1BinnedR,CCGbinsize, Rmaxlags);
globalACG2RawR=CCGraw(globalACG2BinnedR, CCGbinsize, Rmaxlags);

globalACG1SignalR=CCGsignal(globalACG1BinnedR, CCGbinsize, Rmaxlags);
globalACG2SignalR=CCGsignal(globalACG2BinnedR, CCGbinsize, Rmaxlags);

globalACG1NoiseR=globalACG1RawR-globalACG1SignalR;
globalACG2NoiseR=globalACG2RawR-globalACG2SignalR;

globalRRaw=crossCorrelationCoefficient(globalCCGRawR, globalACG1RawR, globalACG2RawR)

globalRSignal=crossCorrelationCoefficient(globalCCGSignalR, globalACG1SignalR, globalACG2SignalR) 

globalRNoise=crossCorrelationCoefficient(globalCCGNoiseR, globalACG1NoiseR, globalACG2NoiseR)


%Plot global PSTH (are plotting PSTH only for cell 1):   
 
figure; plot(0:PSTHbinsize:T, mean(globalBinnedAverageSamples(:,1,:),3)/(PSTHbinsize/1000), '*') ;

title('Global Input PSTH for Sample Cell');
xlabel('time (ms)');
ylabel('Firing rate (Hz)');
axis([0 T 0 100]) ;


%Plot sample global voltage trace from last trial (ie. from trial=trialNumber)

figure; plot(t, globalV);
xlabel('time(ms)');
ylabel('voltage(mV)');
axis([0 4*T -80 -63]);
title(['Voltage trace for superficial cell under global input from a sample trial (trial number ' int2str(numberTrials) ')']); 




 