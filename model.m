%% Bol et al Cancellation of Periodic Input via Feedback (New) 
function [V, tn, weightMatrix] =model(cellNumber, lamda, frequency, stimulation, dt, tfinal, randomNoiseTrain, fixedNoiseTrain, c, e, trainingSession, trialIndex)


%SP Cell parameters

Vthreshold=-65; %mV
Vreset=-68.8; %mV 
tauref=0.7;   %Refractory period of cell (ms). 
cutoffFrequency=500; %Cutoff frequency for noise filtering (Hz) 
tauOmega=4900000; %Time constant for potentiating learning rule (ms) 
eta4=3.6e-3; %LTD learning rule parameter
eta2=1.8e-3; %LTD learning rule parameter
LOmega4=100; %LTD learning rule parameter (ms) 
LOmega2=10; %LTD learning rule parameter (ms) 
T=1000/frequency; %period of sinusoidal input (ms) 

if frequency==0;
    
    T=tfinal;

end 

noPeriods=floor(tfinal/T); 
noSegments=100; %Number of parallel fibres synapsing on each SP cell. 
segmentSize=T/noSegments; %difference in lag between successive parallel fibres. 
gGABA=0.14; %GABA channel conductance (mS/(cm)^2)
gLeak=0.14;  %Leak channel conductance (mS/(cm)^2) 
maxgAMPA=0.024; %Maximum AMPA channel conductance per synapse (mS/(cm)^2) 
kappa=0.21;  %Sinusoidal component amplitude (uA/cm^2) 
capacitance=1; %uF/cm^2
EAMPA=0;  %AMPA channel reversal potential (mV) 
EGABA=-68.8; %GABA channel reversal potential (mV) 
ELeak=Vreset; %Leak channel reversal potential (mV) 
IDC=0.313; %Bias current (uA/cm^2)
noiseSD=0.412; %noise standard deviation 
omegaMax=1;  %maximum value of synaptic weights (also the starting value)


%DAP Parameters 

ADAP=0.6;
gammaDAP=0.2;
taurefDAP=0.1; 
BDAP=2;
betaDAP=0.35;
DDAP=0.1;
alphaDAP=10.9; 
tauDAP=1; 
EDAP=3.5;
tauDAPUnits=7; %ms

%Preallocate matrices

t=0:dt:tfinal; %time vector  
ts=repmat([0:segmentSize:(T-segmentSize)]',1, noPeriods+1)+repmat([0:T:tfinal],noSegments,1);  
%A list of each parallel fibre's most active times. 
tn=zeros(1, tfinal); %this vector keeps track of spike times. 
DAP=zeros(1, length(t));  %Depolarizing after potential
b=zeros(1,length(t));  %DAP parameter b
preRectifiedTerm=zeros(1,length(t)); %feedforward SP cell term
rectifiedTerm=zeros(1,length(t));    %rectified feedforward SP cell term 
V=zeros(1,length(t));   %SP cell voltage
refractoryCounter=repmat(999, 1, length(t));   %if 999 at a given time, then 
%SP cell not in refractory period. If SP cell is in refractory period, this
%variable keeps track of how much of the refractory period has elapsed. 

%Initial variable values

V(1)=Vreset;  %mV 
b(1)=0.01; 
bplusmostRecent=0.1;
tn(1)=-1000004; tn(2)=-1000003; tn(3)=-1000002; tn(4)=-1000001; tn(5)=-1000000;  
%the first five spikes are fictitious spikes that have all occurred long
%before t=0. 
lastTnIndex=5; %keeps track of index of most recent spike in vector tn. 
mostRecentBurstEnd=-1000000; %the most recent burst is set to have occurred before t=0. 
DAP(1)=0; %uA/(cm)^2

%Create a fourth order Butterworth filter and filter both the shared and unique noise trains:
 
nyquistFrequency=(1/2)*(1000/dt); %Hz 
[parameter1,parameter2]=butter(4, cutoffFrequency/nyquistFrequency); 

preNormalizedRandomNoise=filter(parameter1, parameter2, randomNoiseTrain);
preNormalizedFixedNoise=filter(parameter1, parameter2, fixedNoiseTrain);

%Normalize both of these noise trains and set to standard deviation noiseSD:

normRandomNoise=noiseSD*preNormalizedRandomNoise./std(preNormalizedRandomNoise);
normFixedNoise=noiseSD*preNormalizedFixedNoise./std(preNormalizedFixedNoise);

%Create SP cell noise train using both shared and unique noise components: 

noiseSPCell=normFixedNoise*sqrt(c)+normRandomNoise*sqrt(1-c);  %uA/(cm)^2 

%Find AMPAOpeningProbability (which includes opening probability for each AMPA synapse) 

if lamda==1

    AMPAOpeningProbability=AMPAProb(frequency, tfinal, dt, stimulation, noiseSPCell, e, noSegments, trialIndex, cellNumber); 

    else
    
    AMPAOpeningProbability=zeros(1,length(t));

end 

%Create weight matrix. 

if trainingSession ==0;  
    
    % cd('/Users/bensimmonds/Documents/MATLAB/Noise Correlation Reduction Model/Weights'); 
    
    load(['Weights/globalWeightMatrixe' num2str(e) num2str(frequency) 'Hz.mat']);
    
        
else 
  
  weightMatrix=repmat(omegaMax,noSegments, 1); 
  
end 
   

for k=1:length(t)-2; 
        
  %If training weights, update segment weights with potentiation rule for time t(k+1). 

    if trainingSession ==1;
    
        weightMatrix=weightMatrix+dt*(omegaMax-weightMatrix)/tauOmega; 

    end 
            
      
        

         %Find rectified feedforward  term. 
    
         preRectifiedTerm(k)=IDC+noiseSPCell(k)+stimulation*kappa*sin(2*pi*frequency*t(k)/1000);
         rectifiedTerm(k)=preRectifiedTerm(k).*(preRectifiedTerm(k)>0); 
       
        %Next, calculate gAMPA, containing the AMPA conductance for each
        %parallel fibre-SP cell synapse. 

        gAMPA=sum(weightMatrix.*AMPAOpeningProbability(:,k)*maxgAMPA);  %mS/(cm)^2
        
   %Now, will calculate voltage V for time t(k+1). NOTE: lamda is always
   %set to 0 for local input so the feedback is turned off. Under global
   %input lamda is set to 1 and the feedback is turned on. 
       
   V(k+1)=V(k)+(1/capacitance)*dt*(DAP(k)-gLeak*(V(k)-ELeak)+rectifiedTerm(k)+lamda*(-gAMPA*(V(k)-EAMPA)...
           -gGABA*(V(k)-EGABA))); 
      
    
  %Check to see if spike, and if so, reset voltage, set
  %refractory counter to 0 and update spike times vector. 
        
   
       it1=V(k+1)>=Vthreshold; 
       lastTnIndex=it1*(lastTnIndex+1)+(1-it1)*lastTnIndex;
      tn(lastTnIndex)=t(k+1)*it1+tn(lastTnIndex)*(1-it1);  
  
     V(k+1)=it1*Vreset+(1-it1)*V(k+1);
     
     %Set refractory counter to 0 if spike occured. 
     
      refractoryCounter(k+1)=(1-it1)*refractoryCounter(k+1);
      
      
        %If training weights, check for recent bursts, and update weights accordingly. 
        
        if trainingSession==1; 
                      
            it2=it1&&(tn(lastTnIndex-3)>mostRecentBurstEnd); 
            
            it3=it2 &&(tn(lastTnIndex)-tn(lastTnIndex-3) <= 45);  % (BIG BURST - 4 spikes within 45 ms) 
            
            if it3;
           
            mostRecentBurstBeginning=tn(lastTnIndex-3); 
            mostRecentBurstEnd=tn(lastTnIndex); 
            
            difference=ts-mostRecentBurstBeginning;
            preRectifiedUpdateTerms=(1-((difference)./LOmega4).^2); 
            rectifiedUpdateTerms=preRectifiedUpdateTerms.*(preRectifiedUpdateTerms>0) ;
            sumOfRectifiedTerms=sum(rectifiedUpdateTerms,2);
          
            weightMatrix=weightMatrix-eta4*weightMatrix.*sumOfRectifiedTerms;
     
%             
            end 
%             
           
            it4=it1&&tn(lastTnIndex-4)>mostRecentBurstEnd;
            
            it5=it4&&(~it3)&&(tn(lastTnIndex-3)-tn(lastTnIndex-4) <= 15)  ; %(SMALL BURST - 2 spikes within 15 ms)
            
            if it5;
                
      
            mostRecentBurstBeginning=tn(lastTnIndex-4);
            mostRecentBurstEnd=tn(lastTnIndex-3);
           
 
            difference=ts-mostRecentBurstBeginning;
            preRectifiedUpdateTerms=(1-((difference)./LOmega2).^2);
            rectifiedUpdateTerms=preRectifiedUpdateTerms.*(preRectifiedUpdateTerms>0);
            sumOfRectifiedTerms=sum(rectifiedUpdateTerms,2);
 
            weightMatrix=weightMatrix-eta2*weightMatrix.*sumOfRectifiedTerms ;
            
            end 
          
                    
        end   
        
        %find b(k+1), to calculate DAP(k+1)  
        
        spikeDetector=dirac(-tn(1:lastTnIndex)+t(k+1));
       
        it6=sum(spikeDetector)>0;
       
           
        if it6
     
            b(k+1)=b(k)+ADAP+BDAP*(b(k))^2;
            
                 
        else 
     
            
        b(k+1)=b(k)+dt*(-b(k)/(tauDAP*tauDAPUnits));
 
     
        end 
        
        %Update bplusmostRecent to b(k+1) if a spike occured at the previous time t(k) 
        
        it8=abs(tn(lastTnIndex)-t(k))<0.00005; 
        
        if it8 
            
        bplusmostRecent=b(k+1);
        
        
        end 
        
  
        %Compute DAP for time t(k+1).
        
       
        if t(k+1)/tauDAPUnits-tn(lastTnIndex)/tauDAPUnits>taurefDAP && tn(lastTnIndex)/tauDAPUnits-tn(lastTnIndex-1)/tauDAPUnits>dendriteRefractoryPeriod(DDAP, EDAP, bplusmostRecent) ;
       firstTerm= sDAP((t(k+1)/tauDAPUnits-tn(lastTnIndex)/tauDAPUnits), betaDAP*bplusmostRecent);
     secondTerm=sDAP((t(k+1)/tauDAPUnits-tn(lastTnIndex)/tauDAPUnits),gammaDAP);
        DAP(k+1)=alphaDAP*(firstTerm-secondTerm);   
    
       
        
        else 
            DAP(k+1)=0;
        end 
      
     

        
   it9=(refractoryCounter(k+1)<tauref || refractoryCounter(k+1)-tauref<0.00005); 
   
   V(k+1)=Vreset*(it9)+V(k+1)*(1-it9);
   refractoryCounter(k+2)=(refractoryCounter(k+1)+dt)*it9+(1-it9)*999;



end

%Eliminate the five fictitious spikes added at the beginning of the spike
%train.Also get rid of all elements in tn equal to 0. 

tn(1)=[]; tn(1)=[]; tn(1)=[]; tn(1)=[]; tn(1)=[]; 
tn(find(tn==0))=[];



end 


