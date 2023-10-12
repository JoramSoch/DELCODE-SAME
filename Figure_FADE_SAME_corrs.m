% DELCODE-SAME, Figure "FADE SAME corrs"
% _
% original paper: Soch*, Richter* et al., HBM, 2022
% Figure 4: Correlations with independent variables, separated by age group
% 
% DELCODE paper: Soch et al., medRxiv, 2023
% Figure 3: Partial correlations of FADE and SAME scores with other indices of cognitive aging
% 
% written   by Joram Soch <Joram.Soch@DZNE.de>, 12/10/2022, 16:32
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 21/02/2023, 14:48 (DELCODE-SAME)
% finalized by Joram Soch <Joram.Soch@DZNE.de>, 10/10/2023, 16:27


clear
close all
shuffle = false;

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load subject info
addpath(strcat(pwd,'/tools/'));
subj_file = 'subjects/subj_covs.mat';
S = load(subj_file);

% load FADE/SAME scores
FADE_file = 'data/fMRI_scores.xls';
[num, txt, raw] = xlsread(FADE_file);
hdr       = raw(1,:);
data      = raw(2:end,:);
subj_ids  = data(:,1);
num_subj  = size(data,1);
clear num txt

% load behavioral data
bhvr_file = 'subjects/bhvr_data.xls';
[num, txt, raw] = xlsread(bhvr_file);
data_hdr  = raw(1,:);
data_bhvr = raw(2:end,:);
clear num txt


%%% Step 2: process data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign subject diagnosis
subj_grps = {'HC', 'SCD', 'MCI', 'AD', 'AD-rel'};
ind_vars  = {'age', 'age_int', 'education_years', 'employment_years', ...
             'height', 'weight', 'BMI', ...
             'MMSE_total', 'NPT_global' 'PACC5_score', ...
             'Abeta38', 'Abeta40', 'Abeta42', ...
             'total_tau', 'phosphotau181', 'ratio_Abeta42_40'};
num_grps  = numel(subj_grps);
num_vars  = numel(ind_vars);
subj_diag = cell(num_subj,1);
subj_cat  = zeros(num_subj,1);
subj_covs = zeros(num_subj,num_vars);
for i = 1:num_subj
    for j = 1:numel(S.subj_ids)
        if strcmp(subj_ids{i},S.subj_ids{j})
            subj_diag(i)   = S.subj_info(j,strcmp(S.covs,'diagnosis'));
            subj_cat(i)    = find(strcmp(subj_grps,subj_diag{i}));
            for k = 1:num_vars
                subj_covs(i,k) = cell2mat(S.subj_info(j,strcmp(S.covs,ind_vars{k})));
            end;
        end;
    end;
end;
N2 = numel(subj_diag);

% assign behavioral data
bhvr_vars = {'Aprime', 'CorrHitRate'};
bhvr      = NaN(num_subj,numel(bhvr_vars));
for i = 1:num_subj
    for j = 1:size(data_bhvr,1)
        if strcmp(subj_ids{i},data_bhvr{j,1})
            for k = 1:numel(bhvr_vars)
                bhvr(i,k) = data_bhvr{j,strcmp(data_hdr,bhvr_vars{k})};
            end;
        end;
    end;
end;

% shuffle subject diagnosis
if shuffle
    rng(1);
    i2 = randperm(N2);
    subj_diag = subj_diag(i2);
    subj_cat  = subj_cat(i2);
    clear i2
end;


%%% Step 3: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collect scores
scores = {'novelty_FADE', 'novelty_SAME', 'memory_FADE', 'memory_SAME'};
titles = scores;
Y2 = zeros(N2,numel(scores));
for j = 1:numel(scores)
    % scores
    Y2(:,j) = cell2mat(data(:,strcmp(hdr,scores{j})));
    % labels
    titles{j}(strfind(scores{j},'_')) = '-';
end;

% collect covariates
labels    = ind_vars;
labels{2} = 'A_prime';
X2 = NaN(N2,num_vars);
for k = 1:num_vars
    % covariates
    X2(:,k) = subj_covs(:,k);
    if strcmp(labels{k},'A_prime')
        X2(:,k) = bhvr(:,1);
    end;
    % labels
    labels{k}(strfind(labels{k},'_')) = ' ';
end;
labels{2}         =  'A-prime';
labels(end-5:end) = {'Abeta 38', 'Abeta 40', 'Abeta 42', 'total tau', 'phospho-tau 181', 'ratio Abeta 42/40'};

% prepare DELCODE design
Z = zeros(N2,num_grps);
for l = 1:num_grps
    Z(subj_cat==l,l) = 1;
end;

% correlate FADE/SAME scores
rho = NaN(num_vars,numel(scores));
pr  = NaN(num_vars,numel(scores));
for j = 1:numel(scores)
    for k = 1:num_vars
        % extract DELCODE subjects
        X_jk = X2(:,k);
        Y_jk = Y2(:,j);
        % compute partial correlation
        [rho(k,j), pr(k,j)] = partialcorr(X_jk, Y_jk, Z, 'type', 'Pearson', 'rows', 'pairwise', 'tail', 'both');
    end;
end;
clear y_ijk x_ijk


%%% Step 4: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot partial correlations
figure('Name', 'FADE-SAME, Fig. 4', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
p_thr = [0.05, 0.05/(num_vars), 0.05/(num_vars*numel(scores))];
cmap  = [ [repmat([0:0.01:1]',[1 2]), ones(101,1)];
          [ones(100,1), repmat([0.99:-0.01:0]',[1 2])] ];

hold on;
imagesc(rho');
caxis([-0.4, +0.4]);
axis ij equal;
axis([(1-0.5), (num_vars+0.5), (1-0.5), (numel(scores)+0.5)]);
colormap(cmap);
colorbar;
set(gca,'Box','On');
set(gca,'FontSize',12);
set(gca,'XTick',[1:num_vars],'XTickLabel',labels,'XTickLabelRotation',90);
set(gca,'YTick',[1:numel(scores)],'YTickLabel',titles);
xlabel('independent variable', 'FontSize', 12);
ylabel('type of fMRI score', 'FontSize', 12);
for j = 1:numel(scores)
    for k = 1:num_vars
        sig_str = repmat('*',[1 sum(pr(k,j)<p_thr)]);
        text(k, j, sprintf('%0.2f %s', rho(k,j), sig_str), 'FontSize', 12, ...
            'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    end;
end;