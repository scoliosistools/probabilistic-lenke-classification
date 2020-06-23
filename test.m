coronalCobb = randi([0,60],1,3);
coronalBendCobb = randi([0,60],1,3);
sagittalCobb = randi([0,60],1,3);

modifiers = 'A':'C';
lumbarModifier = modifiers(randi(numel(modifiers)));

[CT_prob, LM_prob, STM_prob] = probabilisticLenkeClassification(coronalCobb,coronalBendCobb,sagittalCobb,lumbarModifier);
