function [result]=clustering(S, cls_num, gt)

[C] = SpectralClustering(S,cls_num);
[result] = Clustering8Measure(gt,C);     

%[A nmi avgent] = compute_nmi(gt,C); % NMI
%ACC = Accuracy(C,double(gt)); % ACC
%[fscore,precision,recall] = compute_f(gt,C); %fscore,precision,recall
%[AR,RI,MI,HI]=RandIndex(gt,C);  %ART
%[purity]=CalPurity(gt,C);
end
