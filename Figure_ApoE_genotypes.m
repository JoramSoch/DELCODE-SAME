% DELCODE-SAME, Figure "ApoE genotypes"
% _
% original paper: Soch*, Richter* et al., HBM, 2022
% Figure S3: Analysis of FADE scores as a function of ApoE genotype
% 
% DELCODE paper: Soch et al., medRxiv, 2023
% Figure S4: Comparison of ApoE genotypes to population distribution
% 
% written   by Joram Soch <Joram.Soch@DZNE.de>, 03/11/2022, 13:52
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 20/02/2023, 19:43 (DELCODE-SAME)
% finalized by Joram Soch <Joram.Soch@DZNE.de>, 10/10/2023, 16:51


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
subj_ApoE = cell(num_subj,1);
subj_cat  = zeros(num_subj,1);
for i = 1:num_subj
    for j = 1:numel(S.subj_ids)
        if strcmp(subj_ids{i},S.subj_ids{j})
            subj_diag(i) = S.subj_info(j,strcmp(S.covs,'diagnosis'));
            subj_ApoE(i) = S.subj_info(j,strcmp(S.covs,'ApoE_genotype'));
            subj_cat(i)  = find(strcmp(subj_grps,subj_diag{i}));
        end;
    end;
end;

% shuffle subject covariates
if shuffle
    rng(1);
    i2 = randperm(N2);
    subj_diag = subj_diag(i2);
    subj_ApoE = subj_ApoE(i2);
    subj_cat  = subj_cat(i2);
end;


%%% Step 3a: analyze DELCODE data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create table
ApoE_bin = 1*strcmp(subj_ApoE,'E3/E3') + 2*strcmp(subj_ApoE,'E2/E4') + 2*strcmp(subj_ApoE,'E3/E4') + 2*strcmp(subj_ApoE,'E4/E4');
ApoE_inc = find(ApoE_bin>0);
num_grps = numel(subj_grps);
covs = {'site', 'gender', 'diagnosis'};
Y2A  = Y2(ApoE_inc,:);
T2   = table(Y2A(:,1), Y2A(:,2), Y2A(:,3), Y2A(:,4), subj_diag(ApoE_inc), ApoE_bin(ApoE_inc), ...
            'VariableNames', [scores, covs(3), {'ApoE'}]);

% perform ANOVAs
for j = 1:numel(scores)
    m2(j).model = fitlm(T2, sprintf('%s ~ diagnosis*ApoE', scores{j}));
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
    pc = table2array(m2(i).anova(1:end-1,5));
    for j = 1:numel(row2)
        % store F-/p-value
        R2{j,i} = sprintf('F = %1.2f, %s', F(j), pvalstr(pc(j), p_thr, []));
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
for j = 1:numel(scores)
    for k = 1:num_grps
        for l = 1:2
            if l == 1, subj_ind = find(subj_cat==k &  strcmp(subj_ApoE,'E3/E3')); end;
            if l == 2, subj_ind = find(subj_cat==k & (strcmp(subj_ApoE,'E2/E4') | strcmp(subj_ApoE,'E3/E4') | strcmp(subj_ApoE,'E4/E4'))); end;
            N2A(l,k)      = numel(subj_ind);
            Y_vals{l,k,j} = Y2(subj_ind,j);
            Y_mean(l,k,j) = mean( Y_vals{l,k,j} );
            Y_sem(l,k,j)  = std( Y_vals{l,k,j} )/sqrt(N2A(l,k));
        end;
        [pt(k,j), t(k,j), df, stats] = stattest({Y_vals{1,k,j}, Y_vals{2,k,j}}, 'ttest2');
    end;
end;
clear subj_ind df stats


%%% Step 3b: analyze ApoE genotypes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify genotypes
gens = {'E2/E2', 'E2/E3', 'E2/E4', 'E3/E3', 'E3/E4', 'E4/E4'};
Ns   = [245, 300, 244, 206]';
perc = [0.8, 14.7, 3.6, 58.4, 18.4, 4.0; ...
        0.3, 13.9, 2.0, 55.2, 26.9, 1.6; ...
        0.8, 10.2, 2.0, 63.5, 21.3, 2.0; ...
        0.5, 10.0, 3.9, 63.7, 20.9, 0.5];
% Source: Li et al., Behavioral Genetics, 2019, vol. 49, pp. 455-468.
Na   = round(perc.*repmat(Ns,[1 numel(gens)])./100);
prop = sum(Na,1)./(sum(sum(Na)));
clear Ns Na

% calculte frequencies
y = subj_ApoE; 
g = subj_cat;
N = zeros(num_grps,numel(gens));
for j = 1:num_grps
    for k = 1:numel(gens)
        N(j,k) = sum(strcmp(y(g==j),gens{k}));
    end;
end;
f = N./repmat(sum(N,2),[1 size(N,2)]);
f = [f; prop];

% perform chi^2 GOF tests
chi2  = zeros(1,num_grps);
pc    = zeros(1,num_grps);
df    = zeros(1,num_grps);
bins  = [1:numel(gens)];
for j = 1:num_grps
   [h, pc(j), stats] = chi2gof(bins, 'Ctrs', bins, 'Frequency', N(j,:), ...
                                    'Expected', prop*sum(N(j,:)), 'NParams', 0);
    chi2(j) = stats.chi2stat;
    df(j)   = stats.df;
end;
clear h stats


%%% Step 5a: visualize ApoE genotypes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot ApoE genotypes
figure('Name', 'FADE/SAME, Fig. S3A', 'Color', [1 1 1], 'Position', [25 50 1600 900]);
sps  = [2, 4, 5, 6, 3, 1];
cols = 'bcmgyr';

for j = 1:size(f,1)
    subplot(2,3,sps(j)); hold on;
    for k = 1:numel(gens)
        bar(k, f(j,k), cols(k));
    end;
    axis([(1-1), (numel(gens)+1), 0, 7/10]);
    set(gca,'Box','On');
    set(gca,'XTick',[1:numel(gens)],'XTickLabel',gens);
    if j == size(f,1), legend(gens, 'Location', 'West'); end;
    xlabel('genotype', 'FontSize', 12);
    ylabel('frequency', 'FontSize', 12)
    if j  < size(f,1), title(sprintf('%s (N = %d)', subj_grps{j}, sum(N(j,:))), 'FontSize', 16); end;
    if j == size(f,1), title('population distribution', 'FontSize', 16);      end;
    for k = 1:numel(gens)
        if j < size(f,1)
            text(k, 2/3, num2str(N(j,k)), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        else
            text(k, 2/3, sprintf('%0.3f', prop(k)), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        end;
    end;
    if j < size(f,1)
        text((1-0.5), 1/2, sprintf('chi^2 = %2.2f, %s', chi2(j), pvalstr(pc(j), 0.001, [])), ...
             'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Middle');
    end;
end;