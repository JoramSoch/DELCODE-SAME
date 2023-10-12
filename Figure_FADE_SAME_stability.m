% DELCODE-SAME, Figure "FADE SAME stability"
% _
% original paper: Soch*, Richter* et al., HBM, 2022
% Figure 6: Stability of the FADE scores for older subjects as a function of reference sample
% 
% DELCODE paper: Soch et al., medRxiv, 2023
% Figure S3: Stability of the FADE and SAME scores as a function of reference sample
% 
% written   by Joram Soch <Joram.Soch@DZNE.de>, 07/10/2022, 17:38
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 20/02/2023, 18:18 (DELCODE-SAME)
% finalized by Joram Soch <Joram.Soch@DZNE.de>, 12/10/2023, 01:15


clear
close all
shuffle = false;

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load subject info
addpath(strcat(pwd,'/tools/'));
subj_file = 'subjects/subj_covs.mat';
S = load(subj_file);

% load FADE/SAME scores
FADE_files = {'data/fMRI_scores.xls', 'data/fMRI_scores_yFADE_ref.xls'};
for i = 1:numel(FADE_files)
    [num, txt, raw] = xlsread(FADE_files{i});
    hdr             = raw(1,:);
    data{i}         = raw(2:end,:);
end;
clear num txt
subj_ids = data{1}(:,1);
num_subj = size(data{1},1);


%%% Step 2: process data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign subject diagnosis
subj_grps = {'HC', 'SCD', 'MCI', 'AD', 'AD-rel'};
num_grps  = numel(subj_grps);
subj_diag = cell(num_subj,1);
subj_cat  = zeros(num_subj,1);
for i = 1:num_subj
    for j = 1:numel(S.subj_ids)
        if strcmp(subj_ids{i},S.subj_ids{j})
            subj_diag(i) = S.subj_info(j,strcmp(S.covs,'diagnosis'));
            subj_cat(i)  = find(strcmp(subj_grps,subj_diag{i}));
        end;
    end;
end;
N2 = numel(subj_diag);

% shuffle subject diagnosis
if shuffle
    rng(1);
    i2 = randperm(N2);
    subj_diag = subj_diag(i2);
    subj_cat  = subj_cat(i2);
    clear i2
end;


%%% Step 3: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collect FADE/SAME scores
scores = {'novelty_FADE', 'novelty_SAME', 'memory_FADE', 'memory_SAME'};
titles = scores;
Y2 = zeros(N2,numel(scores),numel(data));
for i = 1:numel(data)
    for j = 1:numel(scores)
        Y2(:,j,i) = cell2mat(data{i}(:,strcmp(hdr,scores{j})));
        titles{j}(strfind(scores{j},'_')) = '-';
    end;
end;
Y = Y2;
g = subj_cat;

% correlate FADE/SAME scores
N = sum(repmat(g,[1 num_grps])==repmat([1:num_grps],[N2 1]));
b = zeros(2,size(Y,2),num_grps);
r = zeros(num_grps,numel(scores));
p = zeros(num_grps,numel(scores));
for i = 1:num_grps
    for j = 1:numel(scores)
        b(:,j,i) = ME_GLM(Y(g==i,j,2), [Y(g==i,j,1), ones(N(i),1)], eye(N(i)));
        [r(i,j), p(i,j)] = corr(Y(g==i,j,1), Y(g==i,j,2));
    end;
end;


%%% Step 4: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot FADE/SAME scores
figure('Name', 'FADE-SAME, Fig. 3', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
cols  = getcols(subj_grps);
p_thr = [0.05, 0.05/(num_grps-1), 0.05/((num_grps-1)*numel(scores))];

for j = 1:numel(scores)
    subplot(2,2,j); hold on;
    Yj = [Y(:,j,1), Y(:,j,2)];
    y_min_max = [min(min(Yj)), max(max(Yj))];
    for i = 1:num_grps
        plot(-10, -10, '.k', 'Color', cols(i,:)./255, 'MarkerSize', 10);
    end;
    plot(y_min_max, y_min_max, '-k', 'LineWidth', 1);
    for i = 1:num_grps
        plot(y_min_max, y_min_max*b(1,j,i)+b(2,j,i), '--k', 'Color', cols(i,:)./255, 'LineWidth', 1);
    end;
    for i = 1:num_grps
        plot(Yj(g==i,1), Yj(g==i,2), '.k', 'Color', cols(i,:)./255, 'MarkerSize', 10);
    end;
    axis([y_min_max, y_min_max]);
    set(gca,'Box','On');
    if j == 4, legend(subj_grps, 'Location', 'SouthEast'); end;
    xlabel('reference map from young AiA', 'FontSize', 12);
    ylabel('reference map from yFADE', 'FontSize', 12);
    title(titles{j}, 'FontSize', 16);
    for i = 1:num_grps
        text(y_min_max(1), y_min_max(2)-(i/16)*diff(y_min_max), ...
             sprintf('   r = %0.2f, %s, y = %0.2f x + %0.2f', r(i,j), pvalstr(p(i,j), 0.001), b(1,j,i), b(2,j,i)), ...
             'Color', cols(i,:)./255, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
    end;
end;