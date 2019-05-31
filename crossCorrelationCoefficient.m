function [ R ] = crossCorrelationCoefficient(CCG, ACG1, ACG2)

R=mean(CCG)./(sqrt(mean(ACG1)).*sqrt(mean(ACG2))); 

end

