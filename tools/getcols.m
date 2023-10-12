function cols = getcols(labs)
% _
% Get DELCODE plotting colors


% get number of labels
num_labs = numel(labs);
cols     = zeros(num_labs,3);

% determine RGB triplets
for i = 1:num_labs
    switch labs{i}
        case 'young'
            cols(i,:) = [  0, 176,  80];
        case 'older'
            cols(i,:) = [  0, 176, 240];
        case 'HC'
            cols(i,:) = [  0, 112, 192];
        case 'SCD'
            cols(i,:) = [255, 192,   0];
        case 'MCI'
            cols(i,:) = [255, 102,   0];
        case 'AD'
            cols(i,:) = [192,   0,   0];
        case 'AD-rel'
            cols(i,:) = [112,  48, 160];
    end;
end;