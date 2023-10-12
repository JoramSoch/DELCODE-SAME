% DELCODE-SAME, Table "FADE SAME"
% _
% original paper: Soch*, Richter* et al., HBM, 2022
% Figure 3: Differences of FADE-classic and FADE-SAME score between age groups
% 
% DELCODE paper: Soch et al., medRxiv, 2023
% Figure 2: FADE and SAME scores as a function of fMRI contrast, score type and diagnostic group
% 
% written   by Joram Soch <Joram.Soch@DZNE.de>, 21/09/2022, 22:28 (V1)
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 14/02/2023, 14:18 (DELCODE-SAME)
% finalized by Joram Soch <Joram.Soch@DZNE.de>, 10/10/2023, 16:05


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

% load FADE/SAME scores (AiA)
FADE_file = 'data/fMRI_scores_AiA_subj.xls';
[num, txt, raw] = xlsread(FADE_file);
data_AiA  = raw(2:end,:);
FADE_inds = 4+[1:4];
age_ind   = 4;
clear num txt


%%% Step 2: process data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign age group
age_grps = {'young', 'older', 'middle-aged'};
age      = cell2mat(data_AiA(:,age_ind));
age_gr   = 1*(age<50) + 2*(age>=60) + 3*(age>=50 & age<60);
N1       = sum(age_gr<3);

% assign subject diagnosis
subj_grps = {'HC', 'SCD', 'MCI', 'AD', 'AD-rel'};
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

% merge AiA and DELCODE
scores = {'novelty_FADE', 'novelty_SAME', 'memory_FADE', 'memory_SAME'};
titles = scores;
Y1 = zeros(N1,numel(scores));
Y2 = zeros(N2,numel(scores));
for j = 1:numel(scores)
    Y1(:,j) = cell2mat(data_AiA(age_gr<3,FADE_inds(j)));
    Y2(:,j) = cell2mat(data(:,strcmp(hdr,scores{j})));
    titles{j}(strfind(scores{j},'_')) = '-';
end;
Y = [Y1; Y2];
g = [age_gr(age_gr<3); subj_cat+2];
subj_grps = [age_grps(1:2), subj_grps];
num_grps  = numel(subj_grps);

% average FADE/SAME scores
N      = sum(repmat(g,[1 num_grps])==repmat([1:num_grps],[N1+N2 1]));
Y_vals = cell(num_grps,numel(scores));
Y_mean = zeros(num_grps,numel(scores));
Y_sem  = zeros(num_grps,numel(scores));
t      = zeros(1,numel(scores));
pt     = zeros(1,numel(scores));
F      = zeros(1,numel(scores));
pF     = zeros(1,numel(scores));
t_HC   = zeros(num_grps,numel(scores));
pt_HC  = zeros(num_grps,numel(scores));
for j = 1:numel(scores)
    for k = 1:num_grps
        Y_vals{k,j} = Y(g==k,j);
        Y_mean(k,j) = mean( Y_vals{k,j} );
        Y_sem(k,j)  = std( Y_vals{k,j} )/sqrt(N(k));
        if strcmp(subj_grps{k},'HC')        % g == 3
            t_HC(k,j)  = NaN;
            pt_HC(k,j) = NaN;
        else                                % g ~= 3
            [pt_HC(k,j), t_HC(k,j), df, stats] = stattest({Y(g==k,j), Y(g==3,j)}, 'ttest2');
        end;
        [pt(j), t(j), df, stats] = stattest({Y(g==1,j), Y(g==2,j)}, 'ttest2');
    end;
    [pF(j), F(j), df, stats] = stattest({Y(g>2,j), g(g>2)}, 'anova1');
end;
clear df stats

% perform more statistical tests
tA  = zeros(3,numel(scores));
pA  = zeros(3,numel(scores));
dfA = zeros(3,numel(scores));
for j = 1:numel(scores)
    for k = 1:size(tA,1)
        % HC/SCD vs. MCI/AD
        if k == 1
            [pA(k,j), tA(k,j), dfA(k,j), stats] = stattest({Y(g==3 | g==4,j), Y(g==5 | g==6 | g==7,j)}, 'ttest2');
        end;
        % SCD vs. MCI
        if k == 2
            [pA(k,j), tA(k,j), dfA(k,j), stats] = stattest({Y(g==4,j), Y(g==5,j)}, 'ttest2');
        end;
        % MCI vs. AD
        if k == 3
            [pA(k,j), tA(k,j), dfA(k,j), stats] = stattest({Y(g==5,j), Y(g==6,j)}, 'ttest2');
        end;
    end;
end;    


%%% Step 4: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot FADE/SAME scores
figure('Name', 'FADE-SAME, Fig. 3', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
cols  = getcols(subj_grps);
p_thr = [0.05, 0.05/(num_grps-1), 0.05/((num_grps-1)*numel(scores))];

for j = 1:numel(scores)
    subplot(2,2,j); hold on;
    for k = 1:num_grps
        bar(k, Y_mean(k,j), 'FaceColor', cols(k,:)./255);
        errorbar(k, Y_mean(k,j), Y_sem(k,j), '.k', 'LineWidth', 2, 'CapSize', 15);
    end;
    plot([1-(1/3), 2+(1/3)], repmat(-(0.8/10)*min(Y_mean(:,j)-Y_sem(:,j)),[1 2]), '-k', 'LineWidth', 1);
    plot([3-(1/3), num_grps+(1/3)], repmat(-(0.8/10)*min(Y_mean(:,j)-Y_sem(:,j)),[1 2]), '-k', 'LineWidth', 1);
    axis([(1-1), (num_grps+1), [(12/10), -(2/10)]*min(Y_mean(:,j)-Y_sem(:,j))]);
    set(gca,'Box','On');
    set(gca,'FontSize',12);
    set(gca,'XTick',[1:num_grps],'XTickLabel',subj_grps);
    xlabel('subject group', 'FontSize', 12);
    ylabel('fMRI score', 'FontSize', 12);
    title(titles{j}, 'FontSize', 16);
    p_str = pvalstr(pt(j), 0.001);
    text(mean([1 2]), -(1/10)*min(Y_mean(:,j)-Y_sem(:,j)), sprintf('t = %0.2f, %s', t(j), p_str), ...
        'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
    p_str = pvalstr(pF(j), 0.001);
    text(mean([3 num_grps]), -(1/10)*min(Y_mean(:,j)-Y_sem(:,j)), sprintf('F = %0.2f, %s', F(j), p_str), ...
        'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
    for k = 1:num_grps
        if isnan(pt_HC(k,j))
            sig_str = '–';
        elseif pt_HC(k,j) > p_thr(1)
            sig_str = 'n.s.';
        else
            sig_str = repmat('*',[1 sum(pt_HC(k,j)<p_thr)]);
        end;
        text(k, (11/10)*min(Y_mean(:,j)-Y_sem(:,j)), sig_str, 'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    end;
    if j == 1
        for k = 1:num_grps
            text(k, (1/10)*min(Y_mean(:,j)-Y_sem(:,j)), sprintf('N = %d', N(k)), ...
                'Color', [1 1 1], 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        end;
    end;
end;

% display test results
for j = 1:numel(scores)
    fprintf('\n-> %s scores:\n', titles{j});
    fprintf('   - HC/SCD/AD-rel vs. MCI/AD: t_%d = %0.2f, %s.\n', dfA(1,j), tA(1,j), pvalstr(pA(1,j), 0.001));
    fprintf('   - SCD vs. MCI: t_%d = %0.2f, %s.\n', dfA(2,j), tA(2,j), pvalstr(pA(2,j), 0.001));
    fprintf('   - MCI vs. AD : t_%d = %0.2f, %s.\n', dfA(3,j), tA(3,j), pvalstr(pA(3,j), 0.001));
end;
fprintf('\n');