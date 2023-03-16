%---------------------------------------------------------------%
%  This is the main script for getting the list of 1's and 0's  %
%  for correct and incorrect trials for each epoch using a      %
%  state-space model -- Wenbo Tang (Jan 05, 2023)               %
%---------------------------------------------------------------%
clc
clear all
close all
% the resulting performance files were also saved in ./Performance_files/
%%
% set the background probability of correct
backprob = 1/3; %1/3 for combined performance; 0.5 for outbound vs. inbound performance
animals = {'JS34'}; % the animal name
%%
%-----------------------------------------------------
%Filter creation
%--------------------------------------------------------
dayfilter = '1';

% runepochfilter = 'isequal($type, ''run'')'; % for all W-track epochs
runepochfilter = 'isequal($environment, ''wtrN'')'; % for novel environment only; use 'wtr2' for AM2
% runepochfilter = 'isequal($environment, ''wtrF'')'; % for novel environment only; use 'wtr1' for AM2
iterator = 'singleepochanal';

f = createfilter('animal',animals, 'days',dayfilter, 'epochs',runepochfilter, 'iterator', iterator);
f = setfilterfunction(f, 'DFAsj_calcproprewarded', {'linpos', 'pos', 'task'});
f = runfilter(f);

% concatenate the days and get the estimated probability of a correct response
for a = 1:length(f)
    behavperform(1).task = 'wtr1';
    for t = 1:length(f(a).output)
        behavperform(t).reward = [];
        behavperform(t).daytrials = zeros(length(f(a).output{t}),2);
        behavperform(t).dayepoch = zeros(length(f(a).output{t}),2);
        trialind = 1;
        for i = 1:length(f(a).output{t})
            b = f(a).output{t}(i);
            nt = length(b.correct) - 1;
            behavperform(t).reward = [behavperform(t).reward; b.correct];
            behavperform(t).daytrials(i,:) = ...
                            [trialind (trialind + nt)];
            behavperform(t).daytime{i} = b.time;
            behavperform(t).dayepoch(i,:) = f(a).epochs{t}(i,:);
            trialind = trialind + nt + 1;
        end
        % get the estimated probability of a correct run
        behavperform(t).probcorrect = getestprobcorrect(...
            behavperform(t).reward, backprob, 0); 
    end
    
    % save results
    save(sprintf('%s%sbehavperform_combine.mat', f(a).animal{2}, ...
                  f(a).animal{3}),'behavperform');

    clear behavperform;
end


