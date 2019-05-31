function [meanCrossCorrelogram] =CCGraw(binnedTrains, binsize, maxlags)   %binsize in ms

%NOTE: spikeTrain input is the train of all spikes, to use to calculate
%mean firing rates of cells. 


binnedCell1=binnedTrains(:,:,1);
binnedCell2=binnedTrains(:,:,2);


trial=binsize*10^-3*size(binnedTrains,2);  %s

%Preallocate correlogram matrix

correlogram=[]; 

for i=1: size(binnedTrains,1);

f1i=mean(binnedCell1(i,:))/(binsize*10^-3);  %Hz
f2i=mean(binnedCell2(i,:))/(binsize*10^-3); %Hz

xcorrVariable=xcorr((binnedCell1(i,:)-binsize*10^-3*f1i),(binnedCell2(i,:)-binsize*10^-3*f2i), maxlags);
nextCorrelogram=xcorrVariable/((binsize*(10^-3))*trial*f1i);

correlogram=[correlogram;nextCorrelogram] ; 

  
end 

meanCrossCorrelogram=mean(correlogram,1); 


end 

