% DELCODE-SAME, Table "FADE SAME dc"
% _
% original paper: Soch*, Richter* et al., HBM, 2022
% Table 3: Within-subject ANOVAs for FADE-classic and FADE-SAME scores
% 
% DELCODE paper: Soch et al., medRxiv, 2023
% Table S3: Effects of diagnosis group and fMRI contrast on single-value scores
% 
% written   by Joram Soch <Joram.Soch@DZNE.de>, 30/09/2022, 17:52
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 20/02/2023, 18:47 (DELCODE-SAME)
% finalized by Joram Soch <Joram.Soch@DZNE.de>, 12/10/2023, 01:00


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

% shuffle subject diagnoses
if shuffle
    rng(1);
    subj_diag = subj_diag(randperm(N2));
end;


%%% Step 3: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create table
T2 = table(Y2(:,1), Y2(:,2), Y2(:,3), Y2(:,4), subj_diag, ...
          'VariableNames', [scores, {'diagnosis'}]);

% perform ANOVAs
for j = 1:2
    Y_ind       = [0,2]+j;
    m2(j).wsf   = table([1 2]', 'VariableNames', {'contrast'});
    m2(j).model = fitrm(T2, sprintf('%s,%s ~ diagnosis', scores{Y_ind(1)}, scores{Y_ind(2)}), 'WithinDesign', m2(j).wsf);
    m2(j).anova = ranova(m2(j).model, 'WithinModel', 'contrast');
end;
clear Y_ind

% generate ANOVA table
p_thr = 0.001;
cols  = {'FADE', 'SAME'};
row2  = m2(1).anova.Properties.RowNames([2,4,5]);
R2    = cell(numel(row2),numel(cols));
for i = 1:numel(cols)
    % extract F-/p-value
    F = table2array(m2(i).anova([2,4,5],4));
    p = table2array(m2(i).anova([2,4,5],5));
    for j = 1:numel(row2)
        % store F-/p-value
        R2{j,i} = sprintf('F = %1.2f, %s', F(j), pvalstr(p(j), p_thr, []));
    end;
end;
clear F p

% calculate effect sizes
dp2 = zeros(1,numel(scores));
psi = zeros(1,numel(scores));
for i = 1:numel(scores)
    % calculate d-prime (HC vs. AD)
    y_HC   = Y2(subj_cat==1,i);
    y_AD   = Y2(subj_cat==4,i);
    s_HC   = var(y_HC);
    s_AD   = var(y_AD);
    s_pool = sqrt( ( (numel(y_HC)-1)*s_HC + (numel(y_AD)-1)*s_AD )/(numel(y_HC)+numel(y_AD)-2) );
    dp2(i) = (mean(y_HC)-mean(y_AD))/s_pool;
    clear y_HC y_AD s_HC s_AD s_pool
    % calculate psi (without AD-rel)
    num_grps = 4;
    n = zeros(1,num_grps);
    m = zeros(1,num_grps);
    s = zeros(1,num_grps);
    for j = 1:num_grps
        n(j) = sum(subj_cat==j);
        m(j) = mean(Y2(subj_cat==j,i));
        s(j) = var(Y2(subj_cat==j,i));
    end;
    gr_avg = mean( Y2(1:num_grps,i) );
    s_pool = sqrt( sum((n-1).*s)/(sum(n)-num_grps) );
    psi(i) = sqrt( 1/(num_grps-1) * sum( ((m - gr_avg)./s_pool).^2 ) );
end;


%%% Step 4: save results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% concatenate ANOVA table
N_str = sprintf('DELCODE (N = %d)', N2);
tab   =      [cell(1,1), cols];
tab   = [tab; row2,      R2];
    
% write ANOVA table
filename = 'Table_FADE_SAME_dc.xls';
xlswrite(filename, tab);
winopen(filename);

% display effect sizes
fprintf('\n-> Effect sizes (HC vs. AD):\n');
for i = 1:numel(scores)
    fprintf('   - %s: d'' = %1.2f, psi = %1.2f.\n', titles{i}, dp2(i), psi(i));
end;
fprintf('\n');