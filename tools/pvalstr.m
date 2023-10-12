function p_str = pvalstr(p, p_thr, sig_thr)
% _
% Create p-value string from p-value


% number of significant digits
num_dig = -floor(log10(p_thr));

% create p-value string
if p < p_thr
    p_str = sprintf(sprintf('p < %%0.%df', num_dig), p_thr);
else
    p_str = sprintf(sprintf('p = %%0.%df', num_dig), p);
end;

% add significance markers
if nargin < 3 || isempty(sig_thr)
    sig_str = '';
else
    sig_str = repmat('*',[1 sum(p<sig_thr)]);
end;
p_str = sprintf('%s%s', p_str, sig_str);