% DELCODE-SAME, Figure "FADE SAME by age"
% _
% original paper: Soch*, Richter* et al., HBM, 2022
% Figure S5: Analysis of FADE scores as a continuous function of age
% 
% DELCODE paper: Soch et al., medRxiv, 2023
% Figure S6: FADE and SAME scores as a continuous function of age (including original study)
% Figure S7: FADE and SAME scores as a continuous function of age (excluding original study)
% 
% written   by Joram Soch <Joram.Soch@DZNE.de>, 09/11/2022, 13:11
% adapted   by Joram Soch <Joram.Soch@DZNE.de>, 20/02/2023, 22:31 (DELCODE-SAME)
% finalized by Joram Soch <Joram.Soch@DZNE.de>, 12/10/2023, 01:40


clear
close all
shuffle = false;
young   = true;
% young = true  -> Figure S6
% young = false -> Figure S7

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

% load FADE/SAME scores (AiA)
FADE_file = 'data/fMRI_scores_AiA_subj.xls';
[num, txt, raw] = xlsread(FADE_file);
data_AiA  = raw(2:end,:);
FADE_inds = 4+[1:4];
age_ind   = 4;
clear num txt


%%% Step 2: process data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign age group
age_grps = {'young', 'middle-aged', 'older'};
age      = cell2mat(data_AiA(:,age_ind));
age_gr   = 1*(age<50) + 2*(age>=50 & age<60) + 3*(age>=60);
N1       = numel(age);

% assign subject diagnosis
subj_grps = {'HC', 'SCD', 'MCI', 'AD', 'AD-rel'};
num_grps  = numel(subj_grps);
subj_diag = cell(num_subj,1);
subj_age  = zeros(num_subj,1);
subj_cat  = zeros(num_subj,1);
for i = 1:num_subj
    for j = 1:numel(S.subj_ids)
        if strcmp(subj_ids{i},S.subj_ids{j})
            subj_diag(i) = S.subj_info(j,strcmp(S.covs,'diagnosis'));
            subj_age(i)  = S.subj_info{j,strcmp(S.covs,'age')};
            subj_cat(i)  = find(strcmp(subj_grps,subj_diag{i}));
        end;
    end;
end;
N2 = numel(subj_diag);

% shuffle subject diagnosis
if shuffle
    rng(1);
    i2a = randperm(N2);
    subj_diag = subj_diag(i2a);
    subj_cat  = subj_cat(i2a);
    i2b = randperm(N2);
    subj_age  = subj_age(i2b);
    clear i2a i2b
end;


%%% Step 3: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collect scores
scores = {'novelty_FADE', 'novelty_SAME', 'memory_FADE', 'memory_SAME'};
titles = scores;
Y1 = zeros(N1,numel(scores));
Y2 = zeros(N2,numel(scores));
for j = 1:numel(scores)
    % scores
    Y1(:,j) = cell2mat(data_AiA(:,FADE_inds(j)));
    Y2(:,j) = cell2mat(data(:,strcmp(hdr,scores{j})));
    % labels
    titles{j}(strfind(scores{j},'_')) = '-';
end;

% collect covariate
x1 = age;
x2 = subj_age;

% merge AiA and DELCODE
Y = [Y1; Y2];
x = [x1; x2];
g = [age_gr; subj_cat+3];
G = max(g);
N = sum(repmat(g,[1 G])==repmat([1:G],[N1+N2 1]));
subj_grps = [age_grps, subj_grps];
num_grps  = numel(subj_grps);

% specify smoothing
if young, xw = 32; else, xw = 16; end;
xp1 = [10:1:88];
xp2 = [58:1:88];
xb  = [(min(x1)-0.5), (35+0.5), 50, (60-0.5), (max(x2)+0.5)];

% smooth AiA scores
mu1 = zeros(numel(xp1),size(Y1,2));
s21 = zeros(numel(xp1),size(Y1,2));
for j = 1:numel(scores)
    for k = 1:numel(xp1)
        yij = Y1(x1>=(xp1(k)-xw/2) & x1<=(xp1(k)+xw/2), j);
        mu1(k,j) = mean(yij);
        s21(k,j) = var(yij);
    end;
end;
clear yij

% smooth DELCODE scores
mu2 = zeros(numel(xp2),size(Y2,2),max(subj_cat));
s22 = zeros(numel(xp2),size(Y2,2),max(subj_cat));
for i = 1:max(subj_cat)
    for j = 1:numel(scores)
        for k = 1:numel(xp2)
            yijk = Y2(subj_cat==i & x2>=(xp2(k)-xw/2) & x2<=(xp2(k)+xw/2), j);
            mu2(k,j,i) = mean(yijk);
            s22(k,j,i) = var(yijk);
        end;
    end;
end;
clear yijk


%%% Step 4: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot smoothing
figure('Name', 'FADE-SAME, Fig. S5', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
cols = getcols(subj_grps);
cols(2,:) = [0, 192, 192];
cols1     = cols(1:3,:);
cols2     = cols(4:8,:);

for j = 1:numel(scores)
    subplot(2,2,j); hold on;
    % formatting
    if j == 3
        if young
            for k = 1:num_grps, plot(-10, -10, '.k', 'Color', cols(k,:)./255, 'MarkerSize', 10); end;
        else
            for k = 4:num_grps, plot(-10, -10, '.k', 'Color', cols(k,:)./255, 'MarkerSize', 10); end;
        end;
    end;  
    if j == 4
        plot(-10, -10, '.k', 'Color', cols(4,:)./255, 'MarkerSize', 10);
        plot(-10, -10, '-k', 'LineWidth', 1);
        plot(-10, -10, '--k', 'LineWidth', 1);
    end;
    plot([min(xp1)-1, max(xp2)+1], [0, 0], '-k');
    % AiA
    if young
        for k = 1:max(age_gr)
            plot(x(g==k), Y(g==k,j), '.b', 'Color', cols1(k,:)./255, 'MarkerSize', 10);
        end;
        plot(xp1, mu1(:,j), '-k', 'LineWidth', 1);
        plot(xp1, mu1(:,j)+s21(:,j), '--k', 'LineWidth', 1);
        plot(xp1, mu1(:,j)-s21(:,j), '--k', 'LineWidth', 1);
    end;
    % DELCODE
    for k = 1:max(subj_cat)
        plot(x(g==k+3), Y(g==k+3,j), '.b', 'Color', cols2(k,:)./255, 'MarkerSize', 10);
        plot(xp2, mu2(:,j,k), '-k', 'Color', cols2(k,:)./255, 'LineWidth', 1);
        if ~young
            plot(xp2, mu2(:,j,k)+s22(:,j,k), '--k', 'Color', cols2(k,:)./255, 'LineWidth', 1);
            plot(xp2, mu2(:,j,k)-s22(:,j,k), '--k', 'Color', cols2(k,:)./255, 'LineWidth', 1);
        end;
    end;
    % formatting
    if young
        for k = 1:numel(xb)
            plot([xb(k), xb(k)], [-100, +100], ':k', 'LineWidth', 1);
        end;
        xlim([min(xp1)-1, max(xp2)+1]);
    else
        xlim([min(xp2), max(xp2)]);
    end;
    ylim([median(Y(:,j))-3*std(Y(:,j)), median(Y(:,j))+3*std(Y(:,j))]);
    set(gca,'Box','On');
    xlabel('age [yrs]', 'FontSize', 12);
    ylabel('fMRI score', 'FontSize', 12);
    title(titles{j}, 'FontSize', 16);
    if j == 3
        if young
            legend(subj_grps(1:end), 'Location', 'NorthWest');
        else
            legend(subj_grps(4:end), 'Location', 'SouthEast');
        end;
    end;
    if j == 4
        legend({'observations', 'smooth mean', 'smooth variance'}, 'Location', 'SouthWest');
    end;
    if j == 2
        if young
            text((18+35)/2, -1.75, 'young subjects', 'FontSize', 10, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
            text((50+59)/2, -1.75, sprintf('middle-\naged'), 'FontSize', 10, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
            text((60+80)/2, -1.75, 'older subjects', 'FontSize', 10, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        end;
    end;
end;