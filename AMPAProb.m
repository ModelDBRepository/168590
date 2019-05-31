function [ permutedAMPAOpeningProbability ] = AMPAProb(frequency, tfinal, dt, stimulation, noiseSPCell, e, noSegments, trialIndex, cellNumber)

pause on

%Parallel fibre parameters

IDC=0.313; %uA/(cm)^2
noiseSD=0.412; 
kappa=0.21; %uA/(cm)^2
Vthreshold=-65; %mV
Vreset=-68.8; %mV 
tauref=0.7; %ms
capacitance=1; %(uF/cm^2)
T=1000/frequency; %ms
gLeak=1/7; %Leak channel conductance (mS/(cm)^2)
ELeak=Vreset; %Leak channel reversal potential (mV) 
cutoffFrequency=500; %Cutoff frequency for noise filtering (Hz) 

if frequency==0;
    
    T=tfinal;

end

blockLength=10000; %total run time broken into blocks to make program more efficient (ms) 
segmentSize=T/noSegments; %difference in lag between successive parallel fibres.  
lag=(0:segmentSize:T-segmentSize)';  %vector containing the lag of each parallel fibre (ms) 

%AMPA Synapse Parameters

Pmax=1;  
taufall=5.26; %ms

%Preallocate vectors

t=0:dt:tfinal;    %ms 
V=repmat(Vreset, T/segmentSize,blockLength/dt+1); %voltage for current time block (mV)
refractoryCounter=repmat(999, T/segmentSize, blockLength/dt+2);
P=zeros(noSegments,blockLength/dt+1); %vector containing AMPA opening probability for each synapse
%for the current time block 
AMPAOpeningProbability=zeros(noSegments, length(t)); %vector containing AMPA opening probability for each
%synapse for total trial time 
tn=zeros(T/segmentSize,blockLength/dt+1);   %list of spike times for current time block 
preRectifiedTerm=zeros(T/segmentSize,1); 
rectifiedTerm=zeros(T/segmentSize,1); 

%Initial variable values

tn(:,1)=repmat(-11000000.00,T/segmentSize,1); 
lastTnIndex=ones(T/segmentSize,1); 

for conductanceIndex=1:(tfinal/blockLength)
    
%Create unfiltered parallel fibre-specific noise and a 4th order Butterworth filter. Each parallel fibre
%has noise that contains a component  identical to the SP cell noise, and a unique component.

%Build Butterworth filter:

nyquistFrequency=(1/2)*(1000/dt); 
preFilteredParallelFibreNoise=normrnd(0,1,blockLength/dt, T/segmentSize); 
[parameter1,parameter2]=butter(4, cutoffFrequency/nyquistFrequency); 

%Filter parallel fibre-specific noise and renormalize to unit variance:

preNormalizedParallelFibreNoise=filter(parameter1,parameter2,preFilteredParallelFibreNoise)';
parallelFibreNoise=noiseSD*(preNormalizedParallelFibreNoise)./repmat(std(preNormalizedParallelFibreNoise,0,2),1,blockLength/dt); 


%Reset vectors at index (1) to value of index (blockLength/dt+1) of
%previous time block (which will be value Vreset if first time block) 

V(:,1)=V(:,blockLength/dt+1); 
V(:,2:blockLength/dt+1)=Vreset;

P(1)=P(blockLength/dt+1);
P(2:blockLength/dt+1)=0; 

%Set vector tn to all zeros, except for tn(:,1), which will have the time of the last
%spike in each parallel fibre. 

tn(:,1)=tn(lastTnIndex); 
tn(:,2:blockLength/dt+1)=zeros(T/segmentSize,blockLength/dt); 

lastTnIndex=ones(T/segmentSize,1); 

refractoryCounter(:,1)=refractoryCounter(:,blockLength/dt+1); 
refractoryCounter(:,2)=refractoryCounter(:,blockLength/dt+2); 
refractoryCounter(:,3:blockLength/dt+2)=999; 


%Iteratively evaluate V at every point on the time vector for this time
%block. 

    for timeIndex=1:blockLength/dt; 
    
%fibreActive is a binary variable indicating whether a parallel fibre is active or not (i.e. has the beginning
%the input arrived at the parallel fibre-SP cell synapse yet) 

fibreActive=(t(timeIndex+(conductanceIndex-1)*(blockLength/dt))>=T | t(timeIndex+(conductanceIndex-1)*(blockLength/dt))>=lag(:,1)); 

%Find noiseSPCell value for each parallel fibre at current time:

noiseSPCellCurrentTime=zeros(noSegments,1); 

        for segmentIndex=1:noSegments;
        
            if fibreActive(segmentIndex,1); 
        
            noiseSPCellCurrentTime(segmentIndex,1)=noiseSPCell(timeIndex+(conductanceIndex-1)*(blockLength/dt)-floor(lag(segmentIndex,1)/dt));
        
        
            end
        
        end

%Calculate rectified feedforward term: 

 
preRectifiedTerm(:,1)=IDC+stimulation*kappa*sin(frequency*2*pi*(t(timeIndex+(conductanceIndex-1)*(blockLength/dt))/1000-lag(:,1)/1000))...
    +sqrt(1-abs(e))*parallelFibreNoise(:,timeIndex)+sign(e)*sqrt(abs(e))*noiseSPCellCurrentTime(:); 
rectifiedTerm(:,1)=preRectifiedTerm(:,1).*(preRectifiedTerm(:,1)>0);

%Calcualte voltage V at next time point (note that V is equal to Vreset if fibre not yet active):

V(:,timeIndex+1)=fibreActive(:,1).*(V(:,timeIndex)+(1/capacitance)*dt*(-gLeak*(V(:,timeIndex)-ELeak)+rectifiedTerm(:,1)))+(1-fibreActive(:,1)).*Vreset; 

 %Check to see if each parallel fibre spikes. If there is a spike, set
 %update pertinent values: 
        
       it1(:,1)=V(:,timeIndex+1)>=Vthreshold; 
       lastTnIndex(:,1)=it1(:,1).*(lastTnIndex(:,1)+1)+(1-it1(:,1)).*lastTnIndex(:,1);
       lastTnInds=sub2ind([T/segmentSize, length(t)], [1:T/segmentSize]', lastTnIndex(:,1)); 
       tn(lastTnInds)=it1(:,1).*t(timeIndex+(conductanceIndex-1)*(blockLength/dt)+1)+(1-it1(:,1)).*tn(lastTnInds) ; 
       refractoryCounter(:,timeIndex+1)=(1-it1(:,1)).*refractoryCounter(:,timeIndex+1);       
       
%Calculate probability of AMPA channels opening at each parallel fibre-SP
%cell synapse: 
 
      
      P(:,timeIndex+1)=it1(:).*(Pmax)+(1-it1(:)).*Ps(taufall,  P(:,timeIndex), repmat(t(timeIndex+(conductanceIndex-1)*...
          (blockLength/dt)+1), noSegments,1)-tn(lastTnInds)); 
      
       AMPAOpeningProbability(:,(conductanceIndex-1)*(blockLength/dt)+timeIndex+1)=P(:,timeIndex+1); 
       
       %Send V of each parallel fibre to Vreset if a spike present.
       
       V(:,timeIndex+1)=it1(:,1).*Vreset+(1-it1(:,1)).*V(:,timeIndex+1);
       
        
 %Maintain V at Vreset if still in refractory period. 
        
         it2(:,1)=refractoryCounter(:,timeIndex+1)<tauref | abs(refractoryCounter(:,timeIndex+1)-tauref)<0.00005; 
         
   V(:,timeIndex+1)=Vreset.*it2(:,1)+V(:,timeIndex+1).*(1-it2(:,1)); 
   
   %Increase refractory counter if still in refractory period.  
   
   refractoryCounter(:,timeIndex+2)=(refractoryCounter(:,timeIndex+1)+dt).*it2(:,1)+(1-it2(:,1))*999;
        
        
    end 

end 


% 
%     if trialIndex==1 && cellNumber==1; 
%           
%              figure; plot(t(1:blockLength/dt), V(1,1:blockLength/dt)); 
%              axis([0 250 -80 -63]);
%              xlabel('Time(ms)');
%              ylabel('Voltage(mV)');
%              title('Sample Individual Granule Cell Voltage Output Cell 1'); 
%              
%              figure; plot(t(1:blockLength/dt), V(10,1:blockLength/dt)); 
%              axis([0 250 -80 -63]);
%              xlabel('Time(ms)');
%              ylabel('Voltage(mV)');
%              title('Sample Individual Granule Cell Voltage Output Cell 10'); 
%                                
%     end 
    
%Now, shift the AMPAOpeningProbability array vertically, so that the
%parallel fibre at the peak of its sinusoid when t=0 is at the top of the
%aray.

lagMaxParallelFibre=T-1/frequency/4*1000;
             
shift=-floor(lagMaxParallelFibre/segmentSize);

permutedAMPAOpeningProbability=circshift(AMPAOpeningProbability, shift);            

end

