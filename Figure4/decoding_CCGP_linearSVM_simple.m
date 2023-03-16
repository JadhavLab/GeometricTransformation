function [accuracy_training,accuracy_test,accuracy_shuffle,p_value,Model] = decoding_CCGP_linearSVM_simple(dat_all,figopt)
%---------------------------------------------------------------%
%  This is the function for cross-condition generalization
%  performance using linear decoder (linear SVMs)               %
%  -- Wenbo Tang (Jul 13, 2022)                                 %
%---------------------------------------------------------------%

% INPUTS:
%
%    dat_all = all data, with struct:
%    dat_all.Features = Features (3 umap dims x N data points)for training trials
%    dat_all.trajlabel = Trial Labels (1 x N data points) for training trials
%    dat_all.Features_test = Features for testing trials
%    dat_all.trajlabel_test = Trial Labels for testing trials
%
%    figopt = 1, plotting the decision boundary defined by the classifer; 0, not plotting
%
% OUTPUTS:
%
%    accuracy_training = prediction accuracy for training trials
%    accuracy_test = prediction accuracy for testing trials
%    accuracy_shuffle = the 1000 accuracy values from trial-label shuffled data
%    p_value = pvalue from the shuffling procedure
%    Model = SVMs model details used for prediction
%%
%------ SVM Training------%
% re-organize data
dat = [dat_all.Features; dat_all.trajlabel]';
dat_test = [dat_all.Features_test; dat_all.trajlabel_test]';

trajlabel =  dat_all.trajlabel';
dat = sortrows(dat, size(dat,2));
randid = randperm(length(trajlabel)); % randomly permute trial labels


%------ SVM training------%
% using the best model to predict (testing)
opt = ['-s 0 -t 0'];
Model = svmtrain(dat(:,end),dat(:,1:end-1), opt); %training
[predict_label,~, ~] = svmpredict(dat(:,end),dat(:,1:end-1), Model,' -q ');

% recalculate the accuracy based on the data distribution (counter unbalaced data)
accuracy_training = 1/2* (length(find(predict_label == 1 & dat(:,end) == 1))./length(find(dat(:,end) == 1)) + ...
    length(find(predict_label == 2 & dat(:,end) == 2))./length(find(dat(:,end) == 2)));

%------ Testing on new trials based on the model trained ------%
[predict_label_test,~, ~] = svmpredict(dat_test(:,end),dat_test(:,1:end-1), Model,' -q ');
% accurcy for predicting the center well during incorrect trials
accuracy_test = 1/2* (length(find(predict_label_test == 1 & dat_test(:,end) == 1))./length(find(dat_test(:,end) == 1)) + ...
    length(find(predict_label_test == 2 & dat_test(:,end) == 2))./length(find(dat_test(:,end) == 2)));
%%
%---plot decision boundary (hyperplane)---%
if figopt
    X = dat(:,1:end-1);
    w = Model.SVs' * Model.sv_coef;
    b = -Model.rho;
    if (Model.Label(1) == -1)
        w = -w; b = -b;
    end
    y_hat = sign(w'*X' + b);
    sv = full(Model.SVs); %support vectors

    % plot support vectors
    hold on
%     plot3(sv(:,1),sv(:,2),sv(:,3),'ko', 'MarkerSize', 10); % not plotting SVs in the current version

    % plot decision boundary
    xgrid=linspace(min(X(:,1)), max(X(:,1)), 30);
    ygrid=linspace(min(X(:,2)), max(X(:,2)), 30);
    [X, Y]=meshgrid(xgrid, ygrid);
    Z=(-b-w(1)*X-w(2)*Y)/w(3);
    
    surf(X, Y, Z,'FaceColor',[190,220,230]/255, 'FaceAlpha',0.5, 'EdgeColor','none')

end
%%
%------ Shuffling to test significance------%
% Trial-label shuffling to get the chance level
run = 1000;% shuffle 1000 times 
% shuffling loop
for i = 1:run
    randid2 = randperm(length(trajlabel)); % shuffle trial labels
    dat_shuffled = [dat_all.Features; dat_all.trajlabel(randid2)]';% shuffle the labels
    
    randid3 = randperm(length(dat_all.trajlabel_test)); % shuffle trial labels
    dat_test_shuffled = [dat_all.Features_test; dat_all.trajlabel_test(randid3)]';% shuffle the labels
    
    % SVM training
    opt = ['-s 0 -t 0'];
    Model_shuf = svmtrain(dat_shuffled(:,end),dat_shuffled(:,1:end-1), opt); %training

    % SVM testing
    %------ Testing on new data based on the model trained ------%
    [predict_label_test,~, ~] = svmpredict(dat_test_shuffled(:,end),dat_test_shuffled(:,1:end-1), Model_shuf,' -q ');

    accuracy_test_shuf = 1/2* (length(find(predict_label_test == 1 & dat_test_shuffled(:,end) == 1))./length(find(dat_test_shuffled(:,end) == 1)) + ...
    length(find(predict_label_test == 2 & dat_test_shuffled(:,end) == 2))./length(find(dat_test_shuffled(:,end) == 2)));

    accuracy_shuffle(i) = accuracy_test_shuf;
    clear accuracy_shuffle_all
end

% calculate p-value based on shuffled data
p_value =  mean(accuracy_shuffle > accuracy_test);
