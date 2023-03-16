%---------------------------------------------------------------%
%  This is the main script for get decoding error of animal's   %
%  current position using rate maps from N' in a trial-by-trial %
%  basis -- Wenbo Tang (Jan 05, 2023)                           %
%---------------------------------------------------------------%
clc
clear all
close all
%%
famnov_dir = '/Volumes/NovelFamiliar/'; % set directory
animalprefix_list = {'AM2','JS17','JS21','ZT2','JS34'};
day_list = 1;
ep_list = [2,4]; % 2 = N; 4 = F;
CTXHP = 1; % 1, PFC; 0, CA1
cellthresh = 1; % minimal number of cells active in a time bin
%%
for theAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{theAnimal};
    % set animal directory
    dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix);
    for day = day_list
        % get trial information
        trajinfo = loaddatastruct(dir, animalprefix, 'trajinfo', day);
        for ep = ep_list
            disp(['Animal: ',animalprefix,' Epoch:',num2str(ep)]) % track decoding process
            %----find trials-----%
            trajbound = trajinfo{day}{ep}.trajbound;
            trajtime = trajinfo{day}{ep}.trajtime;
            if ep == 1
               trajtime = trajtime(end-10+1:end,:);%last 10 trials
            else
               trajtime = trajtime(1:10,:);% first 10 trials
            end
            
            for tr = 1:length(trajtime(:,1))
                current_trajtime = trajtime(tr,:);
                % get decoding info
                if CTXHP
                    [decodinginfo_seg,decodederror_seg] = decoding_position_PFC_novelfamiliar_trial(animalprefix,day,ep,6,current_trajtime,cellthresh);
                else
                    [decodinginfo_seg,decodederror_seg] = decoding_position_CA1_novelfamiliar_trial(animalprefix,day,ep,6,current_trajtime,cellthresh);
                end
                if theAnimal == 1 
                    decodederror{ep}{tr} = decodederror_seg;
                else
                    decodederror{ep}{tr} = [decodederror{ep}{tr};decodederror_seg];
                end
            end
        end
        clear trajinfo
    end
end
%%
% gather statistics
decodederror_stats = [];
for ep = ep_list
    for tr = 1:10
        decodederror_tr = decodederror{ep}{tr};
        errormedian = nanmedian(decodederror_tr);
        errorstd = 1.486*mad(decodederror_tr,1);
        errornum = length(decodederror_tr);
        decodederror_stats = [decodederror_stats;errormedian,errorstd,errornum];
    end
end
%%
%---- plot results ----%
figure,
tr_num = [-10:-1,1:10];
errorbar(tr_num,decodederror_stats(:,1),decodederror_stats(:,2)./sqrt(decodederror_stats(:,3)-1))
hold on
plot([0,0],[0,70],'k--')
ylim([0,70])
xlim([-7,7])
xlabel('#Trial relative to transition')
ylabel('Decoding error using N2 (cm)')
