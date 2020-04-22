function outp = confusion2PerformanceMetrics(C, groups, printFlag)
%Usage:
%  PerformanceMetrics = confusion2PerformanceMetrics(ConfMatrix, GroupNames, PrintFlag)
%
%Description:
%  This function converts the confusion matrix into preformance metrics
%
%Input variables:
%  ConfMatrix: double, the confusion matrix obtained by doing:
%    [ConfMatrix, GroupNames] = confusionmat(OrigScore, PredScore). Here
%    ConfMatrix(i, j) counts observations "i" predicted as "j"
%  GroupNames: same type as OrigScore and PredScore, gives the ordered
%    names of classification groups.
%  PrintFlag: boolean, true turns off displaying output results.
%
%Output variable:
%  PerformanceMetrics: table containing the following classification
%    performance measures for all of the group categories:
%    SEN: sensitivity
%    FPR: false positive rate
%    SPC: specificity
%    ACC: accuracy
%    PPV: positive predictive value
%    NPR: negative positive rate
%  For details see https://en.wikipedia.org/wiki/Confusion_matrix
%
%See also confusionmat
%
%Author: Dmitri Volfson <dmitri.volfson@takeda.com>
%Edited by Tamas Kiss <kiss.t@wigner.hu>

%% Calculate classification performance metrics
res=[];
%metrics
res.colnames={'SEN','FPR','SPC','ACC','PPV','NPV','FDR'};
%groups
res.rownames=groups;
res.dat=zeros(length(res.rownames),length(res.colnames));
inds=1:length(res.rownames);

for ig=1:length(groups)
    TP = C(ig,ig);
    FN = sum(C(ig,inds~=ig));
    FP = sum(C(inds~=ig,ig));
    TN = sum(sum(C(inds~=ig,inds~=ig)));    
    %SEN
    res.dat(ig,1)=TP/(TP+FN);
    %FPR
    res.dat(ig,2)=FP/(FP+TN);
    %SPC
    res.dat(ig,3)=TN/(FP+TN);
    %ACC
    res.dat(ig,4)=(TP+TN)/(TP+FN+FP+TN);
    %PPV
    res.dat(ig,5)=TP/(TP+FP);
    %NPV
    res.dat(ig,6)=TN/(FN+TN);
    %FDR
    res.dat(ig,7)=FP/(FP+TP);
end

%% Create output table
outp = array2table(res.dat);
outp.Properties.VariableNames = res.colnames;
outp.Properties.RowNames = res.rownames;

%% Print table if required
if printFlag
    disp(outp)
end

end