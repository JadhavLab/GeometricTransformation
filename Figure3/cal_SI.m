%---------------------------------------------------------------%
%  This is the script for calculating trajectory selectiviy     %
%  indices for all cells    -- Wenbo Tang (Jan 06, 2023)        %
%---------------------------------------------------------------%
clc
clear all
close all
%%
animalprefix_list = {'AM2','JS17','JS21','ZT2','JS34'}; % all animals
famnov_dir = '/Volumes/NovelFamiliar'; % set directory
day = 1;
eps = 4; % epochs, 2 = N, 4 = F, 6 = N'
HPCTX = 0; % HPCTX = 1, CA1; HPCTX = 0, PFC
%%
% gather selecivity indices
SI_cellinfo = []; % reset matrix

for theAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{theAnimal};
    dir = sprintf('%s/%s_direct/',famnov_dir,animalprefix);
    % interneuron list
    if strcmp(animalprefix,'AM2')
         exclude_list = [6,9];%AM2
    elseif strcmp(animalprefix,'JS21')
         exclude_list = [6,3;6,6;6,7;6,8;7,1;21,1;22,1;25,1];%JS21
    elseif strcmp(animalprefix,'JS17')
        exclude_list = [6,1;7,1];%JS17
    elseif strcmp(animalprefix,'ZT2')
        exclude_list = [31,1;31,5;32,2;27,1;14,10];%ZT2
    elseif strcmp(animalprefix,'JS34')
        exclude_list = [];%JS34
    end
    %%
    % get cell indices, exclude interneurons

    if HPCTX
       [~, cellidx] = matchidx_acrossep(dir, animalprefix, day, exclude_list);% find the common neurons across eps
    else
       [cellidx, ~] = matchidx_acrossep(dir, animalprefix, day, exclude_list);% find the common neurons across eps
    end
    cellnum = length(cellidx(:,1));


    load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
    
    for ep = eps
        disp([animalprefix,' EP-',num2str(ep)]) % display progress
        
         for c = 1:cellnum
              excell_id = cellidx(c,:);
              % calculate SI(IN)
              linfield = linfields{day}{ep}{excell_id(1)}{excell_id(2)};
              tempos1 = linfield{2}(:,1);
              tempfield1 = linfield{2}(:,5);
              tempos2 = linfield{4}(:,1);
              tempfield2 = linfield{4}(:,5);

              tempfield2 = interp1(tempos2,tempfield2,tempos1,'nearest');

              linfield_cell_in = [tempfield1,tempfield2];
              validid = ~isnan(linfield_cell_in);
              validid = validid(:,1) & validid(:,2);
              valididx = find(validid == 1);

              meanrate1 = nanmean(linfield_cell_in(valididx,1));
              meanrate2 = nanmean(linfield_cell_in(valididx,2));

              peakrate1 = max(linfield_cell_in(valididx,1));
              peakrate2 = max(linfield_cell_in(valididx,2));
              peakrate_in =  max(peakrate1,peakrate2);

              Pref_dir_in = sign((sum(linfield_cell_in(valididx,1)-linfield_cell_in(valididx,2))));
              SImean_in = (meanrate1-meanrate2)./(meanrate1+meanrate2);

              % calculate SI(OUT)
              tempos1 = linfield{1}(:,1);
              tempfield1 = linfield{1}(:,5);
              tempos2 = linfield{3}(:,1);
              tempfield2 = linfield{3}(:,5);

              tempfield2 = interp1(tempos2,tempfield2,tempos1,'nearest');

              linfield_cell_in = [tempfield1,tempfield2];
              validid = ~isnan(linfield_cell_in);
              validid = validid(:,1) & validid(:,2);
              valididx = find(validid == 1);

              meanrate1 = nanmean(linfield_cell_in(valididx,1));
              meanrate2 = nanmean(linfield_cell_in(valididx,2));

              peakrate1 = max(linfield_cell_in(valididx,1));
              peakrate2 = max(linfield_cell_in(valididx,2));
              peakrate_out =  max(peakrate1,peakrate2);

              Pref_dir_out = sign((sum(linfield_cell_in(valididx,1)-linfield_cell_in(valididx,2))));
              SImean_out = (meanrate1-meanrate2)./(meanrate1+meanrate2);

              if peakrate_out > 3 && peakrate_in >3 % spatially-tuned cells
                SI_cellinfo = [SI_cellinfo;Pref_dir_in,SImean_in,Pref_dir_out,SImean_out];
              else
                SI_cellinfo = [SI_cellinfo;nan,nan,nan,nan]; % not valid
              end
        end
    end
end
%%
% plot results
figure,
plot(SI_cellinfo(:,2),SI_cellinfo(:,4),'o')
if HPCTX
   title('CA1')
else
   title('PFC')
end
xlabel('SI(IN)')
ylabel('SI(OUT)')
xlim([-1,1])
ylim([-1,1])

