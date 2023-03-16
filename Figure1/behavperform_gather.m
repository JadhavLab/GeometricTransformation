%---------------------------------------------------------------%
%  This is the main script for gathering behavioral performance %
%  for all the animals  -- Wenbo Tang (Jan 05, 2023)            %
%---------------------------------------------------------------%
clc
clear all
close all
%%
% set directory
denovo_dir = '/Volumes/SingledayExp/';
famnov_dir = '/Volumes/NovelFamiliar/';
%%
% the animal list
animal_list = {'AM2','JS17','JS21','ZT2','JS34'};
%%
% gather trial-by-trial performance data
max_trialnum = 50; % plot the first 50 trials
probcorrect_nov_all = [];
probcorrect_denovo_all = [];

for animal = 1:length(animal_list)
    animalprefix = animal_list{animal};
    current_famnov_dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix);
    current_denovo_dir = sprintf('%s/%s_direct/',denovo_dir,animalprefix);
    %------get the behavperf data-----%
    % performance files can also be found in ./Performance_files/
    load(sprintf('%s%sbehavperform_combine.mat',current_famnov_dir,animalprefix)); % get behavior performance for novel sessions
    probcorrect_nov_all = [probcorrect_nov_all,behavperform.probcorrect(1:max_trialnum,1)];
    

    load(sprintf('%s%sbehavperform_combine.mat',current_denovo_dir,animalprefix)); % get behavior performance for de novo sessions
    probcorrect_denovo_all = [probcorrect_denovo_all,behavperform.probcorrect(1:max_trialnum,1)];
end
%%
% gather session-avarged performance data
probcorrect_novavg_all = [];
probcorrect_famavg_all = [];
probcorrect_denovoavg1st_all = [];
probcorrect_denovoavgEnd_all = [];

for animal = 1:length(animal_list)
    animalprefix = animal_list{animal};
    current_famnov_dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix);
    current_denovo_dir = sprintf('%s/%s_direct/',denovo_dir,animalprefix);
    %------get the behavperf data-----%
    % novel sessions
    load(sprintf('%s%sbehavperform_combine.mat',current_famnov_dir,animalprefix)); % get behavior performance for novel sessions
    for sess = 1:length(behavperform.dayepoch(:,1))
        current_trials_range = behavperform.daytrials(sess,:);
        current_probcorrect = behavperform.probcorrect(current_trials_range(1):current_trials_range(2),1);
        current_probcorrect_avg(sess) = mean(current_probcorrect);
    end
    probcorrect_novavg_all = [probcorrect_novavg_all,current_probcorrect_avg'];
   
    % familiar sessions
    load(sprintf('%s%sbehavperform_combine_w1w2.mat',current_famnov_dir,animalprefix)); % get behavior performance for familiar sessions
    current_trials_range = behavperform.daytrials(2,:); %familiar session is the 2nd session
    current_probcorrect = behavperform.probcorrect(current_trials_range(1):current_trials_range(2),1);
    current_probcorrect_avg = mean(current_probcorrect);
    probcorrect_famavg_all = [probcorrect_famavg_all,current_probcorrect_avg];
    
    % de novo sessions
    load(sprintf('%s%sbehavperform_combine.mat',current_denovo_dir,animalprefix)); % get behavior performance for de novo sessions
    current_trials_range = behavperform.daytrials(1,:); %get the 1st de novo session
    current_probcorrect = behavperform.probcorrect(current_trials_range(1):current_trials_range(2),1);
    current_probcorrect_avg = mean(current_probcorrect);
    probcorrect_denovoavg1st_all = [probcorrect_denovoavg1st_all,current_probcorrect_avg];
    current_trials_range = behavperform.daytrials(end,:); %get the last de novo session
    current_probcorrect = behavperform.probcorrect(current_trials_range(1):current_trials_range(2),1);
    current_probcorrect_avg = mean(current_probcorrect);
    probcorrect_denovoavgEnd_all = [probcorrect_denovoavgEnd_all,current_probcorrect_avg];
end

probcorrect_sessionavg_all = [probcorrect_denovoavg1st_all;probcorrect_denovoavgEnd_all;probcorrect_novavg_all(1,:);...
    probcorrect_famavg_all;probcorrect_novavg_all(2,:)];
%%
% plot results
newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54
             0 1 0.8]; % set colors for different animals
figure('position',[1150 500 800 300])
subplot(121)
h = plot((1:max_trialnum)',probcorrect_nov_all,'linewidth',2);
set(h, {'color'}, num2cell(newcolors,2));
hold on
h = plot((1:max_trialnum)',probcorrect_denovo_all,':','linewidth',2);
set(h, {'color'}, num2cell(newcolors,2));
plot([1,max_trialnum],[1/3,1/3],'k--')
ylim([0,1])
xlim([1,50])
xlabel('Trial no.')
ylabel('Prob. correct')


subplot(122)
h = plot(probcorrect_sessionavg_all);
set(h, {'color'}, num2cell(newcolors,2));
hold on
h = plot(probcorrect_sessionavg_all,'o');
set(h, {'color'}, num2cell(newcolors,2));
plot([1,5],[1/3,1/3],'k--')
ylabel('Prob. correct')
xticks(1:5)
xticklabels({'de novo (1st)','de novo (Final)','N','F','N2'})
ylim([0.3,1])

    
