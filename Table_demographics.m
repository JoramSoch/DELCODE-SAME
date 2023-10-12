% DELCODE-SAME, Table "demographics"
% _
% original paper: Soch*, Richter* et al., HBM, 2021
% Table 1: Demographics of young and older subjects
% 
% DELCODE paper: Soch et al., medRxiv, 2023
% Table 1: Demographic information of participant groups
% 
% written   by Joram Soch <Joram.Soch@DZNE.de>, 26/09/2022, 15:05 (V1)
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 31/01/2023, 12:23 (DELCODE-BMS)
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 02/06/2023, 18:00 (DELCODE-SAME)
% finalized by Joram Soch <Joram.Soch@DZNE.de>, 10/10/2023, 15:19


clear
close all

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load subject info
addpath(strcat(pwd,'/tools/'));
subj_file = 'subjects/subj_covs.mat';
S = load(subj_file);

% load subjects of interest
subj_ids = S.subj_ids;
num_subj = numel(subj_ids);
% all subjects, no exclusion

% specify p-value thresholds
p_thr = [0.05, 0.01, 0.001];


%%% Step 2: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get subject information
subj_grps = {'HC', 'SCD', 'MCI', 'AD', 'AD-rel'};
subj_diag = cell(num_subj,1);
subj_cat  = zeros(num_subj,1);
subj_info = cell(num_subj,numel(S.covs));
for i = 1:num_subj
    for j = 1:numel(S.subj_ids)
        if strcmp(subj_ids{i},S.subj_ids{j})
            subj_diag{i}   = S.subj_info{j,strcmp(S.covs,'diagnosis')};
            subj_cat(i)    = find(strcmp(subj_grps,subj_diag{i}));
            subj_info(i,:) = S.subj_info(j,:);
        end;
    end;
end;

% separate subject groups
num_grps  = numel(unique(subj_diag));
subj_covs = cell(1,num_grps);
for h = 1:num_grps
    subj_covs{h} = subj_info(strcmp(subj_diag,subj_grps{h}),:);
end;

% analyze sample size
N_vals = zeros(1,num_grps);
N_strs = cell(1,num_grps);
for h = 1:num_grps
    N_vals(h) = size(subj_covs{h},1);
    N_strs{h} = sprintf('N = %d', N_vals(h));
end;

% analyze age range & mean age
age_vals = zeros(4,num_grps);
age_strs = cell(3,num_grps);
for h = 1:num_grps
    age1 = cell2mat(subj_covs{h}(:,strcmp(S.covs,'age')));
    age2 = cell2mat(subj_covs{h}(:,strcmp(S.covs,'age_int')));
    age_vals(1,h) = min(age2);
    age_vals(2,h) = max(age2);
    age_vals(3,h) = mean(age1);
    age_vals(4,h) = std(age1);
    age_strs{1,h} = sprintf('%d-%d yrs', age_vals(1,h), age_vals(2,h));
    age_strs{2,h} = sprintf('%2.2f ± %0.2f yrs', age_vals(3,h), age_vals(4,h));
    if h == 1
        age_strs{3,h} = '–';
    else
        y1 = age2;
        y2 = cell2mat(subj_covs{1}(:,strcmp(S.covs,'age')));
        [p, t, df, stats] = stattest({y1, y2}, 'ttest2');
        p_str = pvalstr(p, p_thr(end), [p_thr(1), p_thr(1)/(num_grps-1)]);
        age_strs{3,h} = sprintf('t%d = %0.2f, %s', df, t, p_str);
    end;
end;
y = cell2mat(subj_info(:,strcmp(S.covs,'age')));
[p, F, df, stats] = stattest({y, subj_cat}, 'anova1');
p_str = pvalstr(p, p_thr(end));
age_stat_str = sprintf('F%d,%d = %0.2f, %s', df(1), df(2), F, p_str);
clear age1 age2 y1 y2 y

% analyze gender ratio
sex_grps = {'m','f'};
sex_vals = zeros(numel(sex_grps),num_grps);
sex_strs = cell(2,num_grps);
for h = 1:num_grps
    sex = subj_covs{h}(:,strcmp(S.covs,'gender'));
    for k = 1:numel(sex_grps)
        sex_vals(k,h) = sum(strcmp(sex,sex_grps{k}));
    end;
    sex_strs{1,h} = sprintf('%d/%d %s/%s', sex_vals(1,h), sex_vals(2,h), sex_grps{1}, sex_grps{2});
    if h == 1
        sex_strs{2,h} = '–';
    else
        [p, chi2, df, stats] = stattest({sex_vals(:,[1 h])}, 'chi2ind');
        p_str = pvalstr(p, p_thr(end), [p_thr(1), p_thr(1)/(num_grps-1)]);
        sex_strs{2,h} = sprintf('χ²%d = %0.2f, %s', df, chi2, p_str);
    end;
end;
[p, chi2, df, stats] = stattest({sex_vals}, 'chi2ind');
p_str = pvalstr(p, p_thr(end));
sex_stat_str = sprintf('χ²%d = %0.2f, %s', df, chi2, p_str);
clear sex

% analyze acquistion site
site_grps = [2, 5, 8, 10, 13, 14, 17, 18]; % no fMRI subjects: 11, 16
site_vals = zeros(numel(site_grps),num_grps);
site_strs = cell(3,num_grps);
for h = 1:num_grps
    site = cell2mat(subj_covs{h}(:,strcmp(S.covs,'site_ID')));
    for k = 1:numel(site_grps)
        site_vals(k,h) = sum(site==site_grps(k));
        site_strs{1,h} = sprintf('%s%d / ', site_strs{1,h}, site_vals(k,h));
        site_strs{2,h} = sprintf('%sS%d, ', site_strs{2,h}, site_grps(k));
    end;
    site_strs{1,h} = site_strs{1,h}(1:end-3);
    site_strs{2,h} = site_strs{2,h}(1:end-2);
    if h == 1
        site_strs{3,h} = '–';
    else
        [p, chi2, df, stats] = stattest({site_vals(:,[1 h])}, 'chi2ind');
        p_str = pvalstr(p, p_thr(end), [p_thr(1), p_thr(1)/(num_grps-1)]);
        site_strs{3,h} = sprintf('χ²%d = %0.2f, %s', df, chi2, p_str);
    end;
end;
[p, chi2, df, stats] = stattest({site_vals}, 'chi2ind');
p_str = pvalstr(p, p_thr(end));
site_stat_str = sprintf('χ²%d = %0.2f, %s', df, chi2, p_str);
clear site

% analyze ApoE genotype
ApoE_grps = {'E2/E2','E2/E3','E2/E4','E3/E3','E3/E4','E4/E4','n/a'};
ApoE_vals = zeros(numel(ApoE_grps),num_grps);
ApoE_strs = cell(3,num_grps);
for h = 1:num_grps
    ApoE = subj_covs{h}(:,strcmp(S.covs,'ApoE_genotype'));
    for k = 1:numel(ApoE_grps)
        ApoE_vals(k,h) = sum(strcmp(ApoE,ApoE_grps{k}));
    end;
    if ApoE_vals(end,h) ~= 0
        ApoE_str_add = sprintf(' (%d missing)', ApoE_vals(end,h));
    else
        ApoE_str_add = '';
    end;;
    ApoE_strs{1,h} = sprintf('%d / %d / %d / %d / %d / %d%s', ...
                             ApoE_vals(1,h), ApoE_vals(2,h), ApoE_vals(3,h), ApoE_vals(4,h), ApoE_vals(5,h), ApoE_vals(6,h), ApoE_str_add);
    ApoE_strs{2,h} = sprintf('%s, %s, %s, %s, %s, %s', ...
                             ApoE_grps{1}, ApoE_grps{2}, ApoE_grps{3}, ApoE_grps{4}, ApoE_grps{5}, ApoE_grps{6});
    if h == 1
        ApoE_strs{3,h} = '–';
    else
        [p, chi2, df, stats] = stattest({ApoE_vals(1:end-1,[1 h])}, 'chi2ind');
        p_str = pvalstr(p, p_thr(end), [p_thr(1), p_thr(1)/(num_grps-1)]);
        ApoE_strs{3,h} = sprintf('χ²%d = %0.2f, %s', df, chi2, p_str);
    end;
end;
[p, chi2, df, stats] = stattest({ApoE_vals(1:end-1,:)}, 'chi2ind');
p_str = pvalstr(p, p_thr(end));
ApoE_stat_str = sprintf('χ²%d = %0.2f, %s', df, chi2, p_str);
clear ApoE ApoE_str_add

% analyze MMSE total
MMSE_vals = zeros(2,num_grps);
MMSE_strs = cell(2,num_grps);
for h = 1:num_grps
    MMSE = cell2mat(subj_covs{h}(:,strcmp(S.covs,'MMSE_total')));
    MMSE = MMSE(~isnan(MMSE));
    MMSE_vals(1,h) = mean(MMSE);
    MMSE_vals(2,h) = std(MMSE);
    MMSE_strs{1,h} = sprintf('%2.2f ± %0.2f', MMSE_vals(1,h), MMSE_vals(2,h));
    if h == 1
        MMSE_strs{2,h} = '–';
    else
        y1 = MMSE;
        y2 = cell2mat(subj_covs{1}(:,strcmp(S.covs,'MMSE_total')));
        y2 = y2(~isnan(y2));
        [p, z, df, stats] = stattest({y1, y2}, 'mann-whitney');
        p_str = pvalstr(p, p_thr(end), [p_thr(1), p_thr(1)/(num_grps-1)]);
        MMSE_strs{2,h} = sprintf('z = %0.2f, %s', z, p_str);
    end;
end;
y = cell2mat(subj_info(:,strcmp(S.covs,'MMSE_total')));
[p, chi2, df, stats] = stattest({y(~isnan(y)), subj_cat(~isnan(y))}, 'kruskal-wallis');
p_str = pvalstr(p, p_thr(end));
MMSE_stat_str = sprintf('χ²%d = %0.2f, %s', df(1), chi2, p_str);
clear MMSE y y1 y2

% analyze NPT global
NPT_vals = zeros(2,num_grps);
NPT_strs = cell(2,num_grps);
for h = 1:num_grps
    NPT = cell2mat(subj_covs{h}(:,strcmp(S.covs,'NPT_global')));
    NPT = NPT(~isnan(NPT));
    NPT_vals(1,h) = mean(NPT);
    NPT_vals(2,h) = std(NPT);
    NPT_strs{1,h} = sprintf('%0.2f ± %0.2f', NPT_vals(1,h), NPT_vals(2,h));
    if h == 1
        NPT_strs{2,h} = '–';
    else
        y1 = NPT;
        y2 = cell2mat(subj_covs{1}(:,strcmp(S.covs,'NPT_global')));
        [p, t, df, stats] = stattest({y1, y2}, 'ttest2');
        p_str = pvalstr(p, p_thr(end), [p_thr(1), p_thr(1)/(num_grps-1)]);
        NPT_strs{2,h} = sprintf('t%d = %0.2f, %s', df, t, p_str);
    end;
end;
y = cell2mat(subj_info(:,strcmp(S.covs,'NPT_global')));
[p, F, df, stats] = stattest({y(~isnan(y)), subj_cat(~isnan(y))}, 'anova1');
p_str = pvalstr(p, p_thr(end));
NPT_stat_str = sprintf('F%d,%d = %0.2f, %s', df(1), df(2), F, p_str);
clear NPT y y1 y2

% analyze PACC5 score
PACC_vals = zeros(2,num_grps);
PACC_strs = cell(2,num_grps);
for h = 1:num_grps
    PACC = cell2mat(subj_covs{h}(:,strcmp(S.covs,'PACC5_score')));
    PACC = PACC(~isnan(PACC));
    PACC_vals(1,h) = mean(PACC);
    PACC_vals(2,h) = std(PACC);
    PACC_strs{1,h} = sprintf('%0.2f ± %0.2f', PACC_vals(1,h), PACC_vals(2,h));
    if h == 1
        PACC_strs{2,h} = '–';
    else
        y1 = PACC;
        y2 = cell2mat(subj_covs{1}(:,strcmp(S.covs,'PACC5_score')));
        [p, t, df, stats] = stattest({y1, y2}, 'ttest2');
        p_str = pvalstr(p, p_thr(end), [p_thr(1), p_thr(1)/(num_grps-1)]);
        PACC_strs{2,h} = sprintf('t%d = %0.2f, %s', df, t, p_str);
    end;
end;
y = cell2mat(subj_info(:,strcmp(S.covs,'PACC5_score')));
[p, F, df, stats] = stattest({y(~isnan(y)), subj_cat(~isnan(y))}, 'anova1');
p_str = pvalstr(p, p_thr(end));
PACC_stat_str = sprintf('F%d,%d = %0.2f, %s', df(1), df(2), F, p_str);
clear NPT y y1 y2

% analyze Abeta 42/40 ratio
Ab42_vals = zeros(2,num_grps);
Ab42_strs = cell(2,num_grps);
for h = 1:num_grps
    Ab42 = cell2mat(subj_covs{h}(:,strcmp(S.covs,'ratio_Abeta42_40')));
    Ab42 = Ab42(~isnan(Ab42));
    Ab42_vals(1,h) = mean(Ab42);
    Ab42_vals(2,h) = std(Ab42);
    Ab42_strs{1,h} = sprintf('%0.3f ± %0.3f', Ab42_vals(1,h), Ab42_vals(2,h));
    if h == 1
        Ab42_strs{2,h} = '–';
    else
        y1 = Ab42;
        y2 = cell2mat(subj_covs{1}(:,strcmp(S.covs,'ratio_Abeta42_40')));
        y2 = y2(~isnan(y2));
        [p, z, df, stats] = stattest({y1, y2}, 'mann-whitney');
        p_str = pvalstr(p, p_thr(end), [p_thr(1), p_thr(1)/(num_grps-1)]);
        Ab42_strs{2,h} = sprintf('z = %0.2f, %s', z, p_str);
    end;
end;
y = cell2mat(subj_info(:,strcmp(S.covs,'ratio_Abeta42_40')));
[p, chi2, df, stats] = stattest({y(~isnan(y)), subj_cat(~isnan(y))}, 'kruskal-wallis');
p_str = pvalstr(p, p_thr(end));
Ab42_stat_str = sprintf('χ²%d = %0.2f, %s', df(1), chi2, p_str);
clear Ab42 y y1 y2

% analyze Amyloid positivity
Amyl_grps = {'A+','A–'};
Amyl_vals = zeros(numel(Amyl_grps),num_grps);
Amyl_strs = cell(2,num_grps);
for h = 1:num_grps
    Amyl = cell2mat(subj_covs{h}(:,strcmp(S.covs,'ratio_Abeta42_40')));
    Amyl_vals(1,h) = sum(Amyl<=0.08);
    Amyl_vals(2,h) = sum(Amyl> 0.08);
    Amyl_strs{1,h} = sprintf('%d/%d %s/%s', Amyl_vals(1,h), Amyl_vals(2,h), Amyl_grps{1}, Amyl_grps{2});
    if h == 1
        Amyl_strs{2,h} = '–';
    else
        [p, chi2, df, stats] = stattest({Amyl_vals(:,[1 h])}, 'chi2ind');
        p_str = pvalstr(p, p_thr(end), [p_thr(1), p_thr(1)/(num_grps-1)]);
        Amyl_strs{2,h} = sprintf('χ²%d = %0.2f, %s', df, chi2, p_str);
    end;
end;
[p, chi2, df, stats] = stattest({Amyl_vals}, 'chi2ind');
p_str = pvalstr(p, p_thr(end));
Amyl_stat_str = sprintf('χ²%d = %0.2f, %s', df, chi2, p_str);
clear Amyl


%%% Step 3: save results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% concatenate results table
clear p chi2 t F df stats p_str
tab =      [{' '},                  subj_grps,      {'Statistics'}];
tab = [tab; {'sample size'},        N_strs,         {'–'}];
tab = [tab; {'age range'},        age_strs(1,:),  {'–'}];
tab = [tab; {'mean age'},         age_strs(2,:),  age_stat_str];
tab = [tab; {'    test vs. HC'},  age_strs(3,:),  {' '}];
tab = [tab; {'gender ratio'},     sex_strs(1,:),  sex_stat_str];
tab = [tab; {'    test vs. HC'},  sex_strs(2,:),  {' '}];
% tab = [tab; {'acquisition site'}, site_strs(1,:), site_stat_str];
% tab = [tab; {' '},                site_strs(2,:), {' '}];
% tab = [tab; {'    test vs. HC'},  site_strs(3,:), {' '}];
tab = [tab; {'MMSE total'},         MMSE_strs(1,:), MMSE_stat_str];
tab = [tab; {'    test vs. HC'},    MMSE_strs(2,:), {' '}];
tab = [tab; {'NPT global'},         NPT_strs(1,:),  NPT_stat_str];
tab = [tab; {'    test vs. HC'},    NPT_strs(2,:),  {' '}];
tab = [tab; {'PACC5 score'},        PACC_strs(1,:), PACC_stat_str];
tab = [tab; {'    test vs. HC'},    PACC_strs(2,:), {' '}];
% tab = [tab; {'ApoE genotype'},      ApoE_strs(1,:), ApoE_stat_str];
% tab = [tab; {' '},                  ApoE_strs(2,:), {' '}];
% tab = [tab; {'    test vs. HC'},    ApoE_strs(3,:), {' '}];
tab = [tab; {'Aβ 42/40 ratio'},     Ab42_strs(1,:), Ab42_stat_str];
tab = [tab; {'    test vs. HC'},    Ab42_strs(2,:), {' '}];
tab = [tab; {'Amyloid positivity'}, Amyl_strs(1,:), Amyl_stat_str];
tab = [tab; {'    test vs. HC'},    Amyl_strs(2,:), {' '}];
disp(tab);

% save results table
filename = 'Table_demographics.xls';
xlswrite(filename, tab);
winopen(filename);