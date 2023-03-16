function [accuracy,p_value,accuracy_shuf] = decoding_dichotomy_linearSVM(dat_all)
%---------------------------------------------------------------%
%  This is the function for decoding dichotomies using          %
%  linear decoder (linear SVMs)                                 %
%  -- Wenbo Tang (Jul 13, 2022)                                 %
%---------------------------------------------------------------%

% INPUTS:
%
%    dat_all = all data, with struct:
%    dat_all.Features = Features (3 umap dims x N data points)
%    dat_all.trajlabel = Trial Labels (1 x N data points)
%
% OUTPUTS:
%
%    accuracy = prediction accuracy
%    p_value = pvalue from the trial-label shuffling procedure
%    accuracy_shuf = prediction accuracy for shuffled data
%%
%------ SVM Training------%
% re-organize data
dat = [dat_all.Features; dat_all.trajlabel]';

trajlabel =  dat_all.trajlabel';
dat = sortrows(dat, size(dat,2));
randid = randperm(length(trajlabel)); % randomly permute trial labels

% SVM training with a 4-fold cross-validation
Training = dat(randid,1:end-1);
GroupTraining = dat(randid,end);
fold = 4;%4-fold


% C-SVMs with a linear kernel 
opt = ['-s 0 -t 0', ' -v ',num2str(fold),' -q '];%4 fold cross-validation
Model = svmtrain(GroupTraining, Training, opt); %training
accuracy = max(Model);

%%
%------ Trial-label shuffling to test significance------%
%Shuffling to get the chance level
run = 1000;% shuffle 1000 times 
% shuffling loop
for i = 1:run
    randid2 = randperm(length(trajlabel)); % shuffle trial labels
    dat_shuffled = [dat_all.Features; dat_all.trajlabel(randid2)]';% shuffle the labels
    randid = randperm(length(trajlabel));
    
    % SVM training with a 4-fold cross-validation
    Training_shuffle = dat_shuffled(randid,1:end-1);
    GroupTraining_shuffle = dat_shuffled(randid,end);
    % C-SVMs with a linear kernel 
    opt = ['-s 0 -t 0', ' -v ',num2str(fold),' -q '];%4-fold
    Model_shuf = svmtrain(GroupTraining_shuffle, Training_shuffle, opt); %training
    accuracy_test_shuf = max(Model_shuf); 

    accuracy_shuf(i) = accuracy_test_shuf;
end

% calculate p-value based on shuffled data
p_value =  mean(accuracy_shuf > accuracy);
