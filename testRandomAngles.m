coronalCobb = randi([0,60],1,3);
coronalBendCobb = randi([0,60],1,3);
sagittalCobb = randi([0,60],1,3);

[CT_prob, STM_prob] = probabilisticLenkeClassification(coronalCobb,coronalBendCobb,sagittalCobb);
