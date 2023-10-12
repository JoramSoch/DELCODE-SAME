% DELCODE-SAME, Figure "FADE SAME Abitur"
% _
% original paper: Soch*, Richter* et al., HBM, 2022
% Figure S4B: Analysis of FADE scores as a function of educational status
% 
% DELCODE paper: Soch et al., medRxiv, 2023
% Figure S5: FADE and SAME scores by diagnostic group and educational status
% 
% written   by Joram Soch <Joram.Soch@DZNE.de>, 03/11/2022, 13:52
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 21/02/2023, 11:18 (DELCODE-SAME)
% finalized by Joram Soch <Joram.Soch@DZNE.de>, 12/10/2023, 01:25


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
subj_diag = cell(num_subj,1);
subj_Abi  = cell(num_subj,1);
subj_cat  = zeros(num_subj,1);
for i = 1:num_subj
    for j = 1:numel(S.subj_ids)
        if strcmp(subj_ids{i},S.subj_ids{j})
            subj_diag(i) = S.subj_info(j,strcmp(S.covs,'diagnosis'));
            subj_Abi(i)  = S.subj_info(j,strcmp(S.covs,'Abitur'));
            subj_cat(i)  = find(strcmp(subj_grps,subj_diag{i}));
        end;
    end;
end;

% shuffle subject covariates
if shuffle
    rng(1);
    i2 = randperm(N2);
    subj_diag = subj_diag(i2);
    subj_Abi  = subj_Abi(i2);
    subj_cat  = subj_cat(i2);
end;


%%% Step 3b: analyze DELCODE data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create table
Abi_bin  = 1*strcmp(subj_Abi,'yes') + 2*strcmp(subj_Abi,'no');
Abi_inc  = find(Abi_bin>0);
num_grps = numel(subj_grps);
covs = {'site', 'gender', 'diagnosis'};
Y2A  = Y2(Abi_inc,:);
T2   = table(Y2A(:,1), Y2A(:,2), Y2A(:,3), Y2A(:,4), subj_diag(Abi_inc), Abi_bin(Abi_inc), ...
            'VariableNames', [scores, covs(3), {'Abitur'}]);

% perform ANOVAs
for j = 1:numel(scores)
    m2(j).model = fitlm(T2, sprintf('%s ~ diagnosis*Abitur', scores{j}));
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
N2A    = zeros(2,num_grps);
Y_vals = cell(2,num_grps,numel(scores));
Y_mean = zeros(2,num_grps,numel(scores));
Y_sem  = zeros(2,num_grps,numel(scores));
t      = zeros(num_grps,numel(scores));
pt     = zeros(num_grps,numel(scores));
df     = zeros(num_grps,numel(scores));
for j = 1:numel(scores)
    for k = 1:num_grps
        for l = 1:2
            if l == 1, subj_ind = find(subj_cat==k & strcmp(subj_Abi,'yes')); end;
            if l == 2, subj_ind = find(subj_cat==k & strcmp(subj_Abi,'no'));  end;
            N2A(l,k)      = numel(subj_ind);
            Y_vals{l,k,j} = Y2(subj_ind,j);
            Y_mean(l,k,j) = mean( Y_vals{l,k,j} );
            Y_sem(l,k,j)  = std( Y_vals{l,k,j} )/sqrt(N2A(l,k));
        end;
        [pt(k,j), t(k,j), df(k,j), stats] = stattest({Y_vals{1,k,j}, Y_vals{2,k,j}}, 'ttest2');
    end;
end;
clear subj_ind stats


%%% Step 4a: save ANOVA results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % concatenate results table
% N_str = sprintf('DELCODE (N = %d)', size(Y2,1));
% tab   =      [cell(1,1), cols];
% tab   = [tab; row2,      R2];
% 
% % save results table
% filename = 'Figure_FADE_SAME_Abitur.xls';
% xlswrite(filename, tab);
% winopen(filename);


%%% Step 4b: visualize Abitur effects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot FADE/SAME scores
figure('Name', 'FADE/SAME, Fig. S4B', 'Color', [1 1 1], 'Position', [75 50 1600 900]);
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
    text(mean([1 (2*num_grps)]), (-2+0.1), sprintf('Abitur main effect: %s', R2{2,j}), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    for k = 1:num_grps
        if pt(k,j) > 0.05, sig_str = 'n.s.';
        else, sig_str = repmat('*',[1 sum(pt(k,j)<p_thr)]); end;
        text((k-1)*2+1.5, 0.1, sig_str, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    end;
    if j == 1
        for k = 1:num_grps
            for l = 1:2
                text((k-1)*2+l, (0-0.1), sprintf('%d', N2A(l,k)), 'Color', [1 1 1], 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
            end;
        end;
    end;    
    if j == 2
        legend(repmat({'with Abitur', 'w/o Abitur'}, [1 num_grps]), 'Location', 'SouthWest');
    end;
end;