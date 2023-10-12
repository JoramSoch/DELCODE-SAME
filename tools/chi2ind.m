function [h, p, cc, stats] = chi2ind(Y, corr);
% _
% Chi-square test for independence given count data
% 
% Sources:
% - https://de.wikipedia.org/wiki/Chi-Quadrat-Test#Unabh%C3%A4ngigkeitstest
% - https://de.wikipedia.org/wiki/Kontingenzkoeffizient#Kontingenzkoeffizient_nach_Karl_Pearson


% get observed counts
O = Y;

% get expected counts
N = sum(sum(O));
E = (sum(O,2)*sum(O,1))/N;

% assess statistical significance
alpha      = 0.05;
stats.chi2 = sum(sum( ((O-E).^2)./E ));
stats.df   = (size(Y,1)-1)*(size(Y,2)-1);
p          = 1-chi2cdf(stats.chi2, stats.df);
h          = (p<alpha);

% calculate contingency coefficient
cc = sqrt(stats.chi2/(stats.chi2 + N));
if strcmp(corr,'yes')
    k  = min([size(Y,1), size(Y,2)]);
    cc = sqrt(k/(k-1)) * cc;
end;