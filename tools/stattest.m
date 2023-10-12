function [p, tstat, df, stats] = stattest(data, test, settings)
% _
% Perform statistical test for measured data


%%% chi-square test for indepence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test,'chi2ind')
    Y = data{1};
    if nargin < 3 || isempty(settings)
        settings.corr = 'yes';
    else
        if ~isfield(settings,'corr'), settings.corr = 'yes'; end;
    end;
    [h, p, cc, stats] = chi2ind(Y, settings.corr);
    tstat = stats.chi2;
    df    = stats.df;
    stats.cc = cc;
end;


%%% two-sample t-test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test,'ttest2')
    y1 = data{1};
    y2 = data{2};
    if nargin < 3 || isempty(settings)
        settings.tail    = 'both';
        settings.vartype = 'equal';
    else
        if ~isfield(settings,'tail'),    settings.tail    = 'both';  end;
        if ~isfield(settings,'vartype'), settings.vartype = 'equal'; end;
    end;
    [h, p, ci, stats] = ttest2(y1, y2, 'tail', settings.tail, 'vartype', settings.vartype);
    tstat = stats.tstat;
    df    = stats.df;
    stats.ci = ci;
end;


%%% one-way ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test,'anova1')
    y = data{1};
    g = data{2};
    if nargin < 3 || isempty(settings)
        settings.display = 'off';
    else
        if ~isfield(settings,'display'), settings.display = 'off'; end;
    end;
    [p, tab, stats] = anova1(y, g, settings.display);
    tstat = tab{2,5};
    df    = [tab{2,3}, tab{3,3}];
end;


%%% rank-sum test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test,'mann-whitney')
    y1 = data{1};
    y2 = data{2};
    if nargin < 3 || isempty(settings)
        settings.method = 'default';
        settings.tail   = 'both';
    else
        if ~isfield(settings,'method'), settings.method = 'default'; end;
        if ~isfield(settings,'tail'),   settings.tail   = 'both';    end;
    end;
    if strcmp(settings.method,'default')
        [p, h, stats] = ranksum(y1, y2, 'tail', settings.tail);
    else
        [p, h, stats] = ranksum(y1, y2, 'method', settings.method, 'tail', settings.tail);
    end;
    tstat = stats.zval;
    df    = NaN;
end;


%%% one-way ANOVA on ranks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test,'kruskal-wallis')
    y = data{1};
    g = data{2};
    if nargin < 3 || isempty(settings)
        settings.display = 'off';
    else
        if ~isfield(settings,'display'), settings.display = 'off'; end;
    end;
    [p, tab, stats] = kruskalwallis(y, g, settings.display);
    tstat = tab{2,5};
    df    = [tab{2,3}, tab{3,3}];
end;