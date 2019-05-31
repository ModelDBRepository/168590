function [ meanCrossCorrelogram] = CCGsignal( binnedTrains,  binsize, maxlags)   %binsize in ms

binnedCell1=binnedTrains(:,:,1);
binnedCell2=binnedTrains(:,:,2);

trial=binsize*10^-3*size(binnedTrains,2); %s

correlogram=[]; 


for i=1:size(binnedTrains,1);

f1i=mean(binnedCell1(i,:))/(binsize*10^-3); %Hz

    for j=1:i-1;
        
        f2j=mean(binnedCell2(j,:))/(binsize*10^-3); %Hz 

        nextCorrelogram=xcorr(binnedCell1(i,:),binnedCell2(j,:), maxlags)/(binsize*(10^-3)*trial*f1i)-f2j; 
     
        correlogram=[correlogram; nextCorrelogram]; 
        
    end 
end 
    
meanCrossCorrelogram=mean(correlogram,1); 


end 



