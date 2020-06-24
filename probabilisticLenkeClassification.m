function [CT_prob, STM_prob] = probabilisticLenkeClassification(coronalCobb,coronalBendCobb,sagittalCobb,CoronalCobb_SD,SagittalCobb_SD)
%PROBABILISTICLENKECLASSIFICATION: calculates the probability of each Lenke
%classification, given the inherent interobserver variability
%
% estimating the probability of each curve type / sagittal thoracic modifier 
% assuming the measured Cobb lies in the center of normal distribution which 
% represents the true underlying Cobb
%
% inputs:
% coronalCobb(1) = PT angle
% coronalCobb(2) = MT angle
% coronalCobb(3) = TL/L angle
% coronalBendCobb(1) = PT bend angle
% coronalBendCobb(2) = MT bend angle
% coronalBendCobb(3) = TL/L bend angle
% sagittalCobb(1) = T2-T5 kyphosis angle
% sagittalCobb(2) = T5-T12 kyphosis angle
% sagittalCobb(3) = T10-L2 kyphosis angle
%
% standard deviation of estimated Cobb angle around the 'true' underlying 
% angle divided into different SD for coronal and sagittal, and different
% SD for different angles:
% CoronalCobb_SD = standard deviation in interobserver variability of coronal cobb angles
% SagittalCobb_SD = standard deviation in interobserver variability of sagittal cobb angles
% Could be further divided into SD for each individual angle - e.g. coronal
% bend PT SD = ...
%
% outputs:
% CT_prob = probability of each curve-type (1-6 -> 1-6)
% STM_prob = probability of each sagittal thoracic modifier (1,2,3 -> -,n,+)

% initialise outputs
CT_prob = zeros(1, 6); % curve-type probabilities
LM_prob = zeros(1, 3); % lumbar modifier probabilities
STM_prob = zeros(1, 3); % sagittal thoracic modifier probabilities

%% Curve Type
p_tl_major = probAGreaterThanB(coronalCobb(3),coronalCobb(2),CoronalCobb_SD)*probAGreaterThanB(coronalCobb(3),coronalCobb(1),CoronalCobb_SD);
p_mt_major = 1 - p_tl_major;

% coronal angle probabilities
p_coronalLess25 = zeros(1, 3);
count = 1;
for ang = coronalCobb
    p_coronalLess25(count) = normcdf(25,ang,CoronalCobb_SD);
    count = count + 1;
end
p_coronalGreater25 = 1 - p_coronalLess25;

% coronal bend angle probabilities
p_coronalBendLess25 = zeros(1, 3);
count = 1;
for ang = coronalBendCobb
    p_coronalBendLess25(count) = normcdf(25,ang,CoronalCobb_SD);
    count = count + 1;
end
p_coronalBendGreater25 = 1 - p_coronalBendLess25;

% sagittal angle probabilities
p_sagittalLess20 = zeros(1, 3);
count = 1;
for ang = sagittalCobb
    p_sagittalLess20(count) = normcdf(20,ang,SagittalCobb_SD);
    count = count + 1;
end
p_sagittalGreater20 = 1 - p_sagittalLess20;

% minor angle is nonstructural if all angles are less than the threshold
% OR if only the coronal (not bend) angle is greater than the threshold
p_minorStructural = (1 - ((p_coronalLess25.* p_coronalBendLess25.*p_sagittalLess20) + (p_coronalGreater25.* p_coronalBendLess25.*p_sagittalLess20)));
p_minorNonStructural = 1 - p_minorStructural;

CT_prob(1) = p_minorNonStructural(1) * p_mt_major * p_minorNonStructural(3);
CT_prob(2) = p_minorStructural(1) * p_mt_major * p_minorNonStructural(3);
CT_prob(3) = p_minorNonStructural(1) * p_mt_major * p_minorStructural(3);
CT_prob(4) = (p_minorStructural(1) * p_mt_major * p_minorStructural(3)) + (p_minorStructural(1) * p_minorStructural(2) * p_tl_major);
CT_prob(5) = p_minorNonStructural(1) * p_minorNonStructural(2) * p_tl_major;
CT_prob(6) = p_minorNonStructural(1) * p_minorStructural(2) * p_tl_major;

%% Sagittal thoracic modifier
STM_prob(1) = normcdf(10,sagittalCobb(2),SagittalCobb_SD);
STM_prob(3) = 1 - normcdf(40,sagittalCobb(2),SagittalCobb_SD);
STM_prob(2) = 1 - (STM_prob(1) + STM_prob(3));

end


function [prob] = probAGreaterThanB(angA, angB, cobb_SD)
    if cobb_SD ~= 0
        step = 0:0.01:180;
        angA_pdf = normpdf(step,angA,cobb_SD);
        angB_pdf = normpdf(step,angB,cobb_SD);
        if angA > angB
            prob = (1-(trapz(min(angA_pdf, angB_pdf))/(2*trapz(angA_pdf))));
        else
            prob = (trapz(min(angA_pdf, angB_pdf))/(2*trapz(angA_pdf)));
        end
    else
        if angA > angB
            prob = 1;
        elseif angA == angB
            prob = 0.5;
        else
            prob = 0;
        end
    end
end

