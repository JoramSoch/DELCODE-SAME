% DELCODE-SAME, Table "FADE SAME sgd"
% _
% original paper: Soch*, Richter* et al., HBM, 2022
% Table 2: Between-subject ANOVAs for FADE-classic and FADE-SAME scores
% 
% DELCODE paper: Soch et al., medRxiv, 2023
% Figure 2: Effects of site, gender and diagnosis on fMRI scores
% 
% written   by Joram Soch <Joram.Soch@DZNE.de>, 30/09/2022, 17:16
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 14/02/2023, 15:55 (DELCODE-SAME)
% finalized by Joram Soch <Joram.Soch@DZNE.de>, 10/10/2023, 16:15


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
hdr  = raw(1,:);
data = raw(2:end,:);
subj_ids = data(:,1);
num_subj = size(data,1);
clear num txt


%%% Step 2: process data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract FADE/SAME scores
scores = {'novelty_FADE', 'novelty_SAME', 'memory_FADE', 'memory_SAME'};
titles = scores;
Y2 = zeros(num_subj,numel(scores));
for j = 1:numel(scores)
    Y2(:,j) = cell2mat(data(:,strcmp(hdr,scores{j})));
    titles{j}(strfind(scores{j},'_')) = '-';
end;
N2 = size(Y2,1);

% assign subject covariates
subj_grps = {'HC', 'SCD', 'MCI', 'AD', 'AD-rel'};
subj_site = zeros(num_subj,1);
subj_sex  = cell(num_subj,1);
subj_diag = cell(num_subj,1);
for i = 1:num_subj
    for j = 1:numel(S.subj_ids)
        if strcmp(subj_ids{i},S.subj_ids{j})
            subj_site(i) = S.subj_info{j,strcmp(S.covs,'site_ID')};
            subj_sex(i)  = S.subj_info(j,strcmp(S.covs,'gender'));
            subj_diag(i) = S.subj_info(j,strcmp(S.covs,'diagnosis'));
        end;
    end;
end;
subj_site = cellstr(num2str(subj_site));

% shuffle subject covariates
if shuffle
    rng(1);
    subj_site = subj_site(randperm(N2));
    subj_sex  = subj_sex(randperm(N2));
    subj_diag = subj_diag(randperm(N2));
end;


%%% Step 3: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create table
covs = {'site', 'gender', 'diagnosis'};
T2   = table(Y2(:,1), Y2(:,2), Y2(:,3), Y2(:,4), subj_site, subj_sex, subj_diag, ...
            'VariableNames', [scores, covs]);

% perform ANOVAs
for j = 1:numel(scores)
    m2(j).model = fitlm(T2, sprintf('%s ~ site + gender*diagnosis', scores{j}));
    m2(j).anova = anova(m2(j).model);
end;

% generate results table
p_thr = 0.001;
cols  = titles;
row2  = m2(1).anova.Properties.RowNames(1:end-1);
R2    = cell(numel(row2),numel(cols));
for i = 1:numel(cols)
    % extract F-/p-value
    F = table2array(m2(i).anova(1:end-1,4));
    p = table2array(m2(i).anova(1:end-1,5));
    for j = 1:numel(row2)
        % store F-/p-value
        R2{j,i} = sprintf('F = %1.2f, %s', F(j), pvalstr(p(j), p_thr, []));
    end;
end;
clear F p


%%% Step 4: save results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% concatenate results table
N_str = sprintf('DELCODE (N = %d)', N2);
tab   =      [cell(1,1), cols];
tab   = [tab; row2,      R2];
    
% save results table
filename = 'Table_FADE_SAME_sgd.xls';
xlswrite(filename, tab);
winopen(filename);