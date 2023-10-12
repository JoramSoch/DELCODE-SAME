% DELCODE-SAME, Figure "FADE SAME gender"
% _
% original paper: Soch*, Richter* et al., HBM, 2022
% Figure: not available
% 
% DELCODE paper: Soch et al., medRxiv, 2023
% Figure S2: FADE and SAME scores by diagnostic group and participant gender
% 
% written   by Joram Soch <Joram.Soch@DZNE.de>, 03/11/2022, 13:52
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 20/02/2023, 14:16 (DELCODE-SAME)
% finalized by Joram Soch <Joram.Soch@DZNE.de>, 12/10/2023, 00:50


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
num_grps  = numel(subj_grps);
subj_diag = cell(num_subj,1);
subj_sex  = cell(num_subj,1);
subj_cat  = zeros(num_subj,1);
for i = 1:num_subj
    for j = 1:numel(S.subj_ids)
        if strcmp(subj_ids{i},S.subj_ids{j})
            subj_diag(i) = S.subj_info(j,strcmp(S.covs,'diagnosis'));
            subj_sex(i)  = S.subj_info(j,strcmp(S.covs,'gender'));
            subj_cat(i)  = find(strcmp(subj_grps,subj_diag{i}));
        end;
    end;
end;

% shuffle subject covariates
if shuffle
    rng(1);
    i2 = randperm(N2);
    subj_diag = subj_diag(i2);
    subj_sex  = subj_sex(i2);
    subj_cat  = subj_cat(i2);
end;


%%% Step 3: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create table
covs = {'site', 'gender', 'diagnosis'};
T2   = table(Y2(:,1), Y2(:,2), Y2(:,3), Y2(:,4), subj_sex, subj_diag, ...
            'VariableNames', [scores, covs([2,3])]);

% perform ANOVAs
for j = 1:numel(scores)
    m2(j).model = fitlm(T2, sprintf('%s ~ diagnosis*gender', scores{j}));
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

% average FADE/SAME scores
N2G    = zeros(2,num_grps);
Y_vals = cell(2,num_grps,numel(scores));
Y_mean = zeros(2,num_grps,numel(scores));
Y_sem  = zeros(2,num_grps,numel(scores));
t      = zeros(num_grps,numel(scores));
pt     = zeros(num_grps,numel(scores));
for j = 1:numel(scores)
    for k = 1:num_grps
        for l = 1:2
            if l == 1, subj_ind = find(subj_cat==k & strcmp(subj_sex,'m')); end;
            if l == 2, subj_ind = find(subj_cat==k & strcmp(subj_sex,'f')); end;
            N2G(l,k)      = numel(subj_ind);
            Y_vals{l,k,j} = Y2(subj_ind,j);
            Y_mean(l,k,j) = mean( Y_vals{l,k,j} );
            Y_sem(l,k,j)  = std( Y_vals{l,k,j} )/sqrt(N2G(l,k));
        end;
        [pt(k,j), t(k,j), df, stats] = stattest({Y_vals{1,k,j}, Y_vals{2,k,j}}, 'ttest2');
    end;
end;
clear subj_ind df stats


%%% Step 4: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot FADE/SAME scores
figure('Name', 'FADE/SAME, Tab. 2', 'Color', [1 1 1], 'Position', [75 50 1600 900]);
cols  = getcols(subj_grps);
p_thr = [0.05, 0.05/num_grps, 0.05/(num_grps*numel(scores))];

for j = 1:numel(scores)
    subplot(2,2,j); hold on;
    for k = 1:num_grps
        for l = 1:2
            bar((k-1)*2+l, Y_mean(l,k,j), 'FaceColor', (3/(2+l))*cols(k,:)./255);
        end;
    end;
    errorbar([1:(2*num_grps)], reshape(Y_mean(:,:,j),[1 (2*num_grps)]), reshape(Y_sem(:,:,j),[1 (2*num_grps)]), ...
             '.k', 'LineWidth', 2, 'CapSize', 15);
    xlim([(1-0.5), (2*num_grps+0.5)]);
    ylim([-2, +0.2]);
    set(gca,'Box','On');
    set(gca,'XTick',[(1+0.5):2:(2*num_grps)],'XTickLabel',subj_grps,'XTickLabelRotation',0);
    xlabel('subject group', 'FontSize', 12);
    ylabel('fMRI score', 'FontSize', 12);
    title(titles{j}, 'FontSize', 16);
  % text(mean([1 (2*num_grps)]), (-2+0.1), sprintf('gender main effect: %s', R2{1,j}), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    for k = 1:num_grps
        if pt(k,j) > 0.05, sig_str = 'n.s.';
        else, sig_str = repmat('*',[1 sum(pt(k,j)<p_thr)]); end;
        text((k-1)*2+1.5, 0.1, sig_str, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    end;
    if j == 1
        for k = 1:num_grps
            for l = 1:2
                text((k-1)*2+l, (0-0.1), sprintf('%d', N2G(l,k)), 'Color', [1 1 1], 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
            end;
        end;
    end;    
    if j == 2
        legend(repmat({'male', 'female'}, [1 num_grps]), 'Location', 'SouthEast');
    end;
end;