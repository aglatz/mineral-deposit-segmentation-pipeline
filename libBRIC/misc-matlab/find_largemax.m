% Takes the output of idx=find(x == max(x)), where x is a vector, and
% returns an idx within the widest maximum as specified by pos.
% INPUTS: idx - value returned by 'find(x == max(x))', where x is a vector.
%         pos - Determines the index that should be returned. Possible
%               values are 'start', 'end', '' or empty, which indicate that
%               the first, last or middle index of the largest maximum
%               should be returned.
% OUTPUTS: idx_max - an index within the largest maximum of x
%          len - width of the largest maximum
% EXAMPLE:
%          idx = [2 3 8 9 10 33];
%          find_largemax(idx, '') % returns 9
%          find_largemax(idx, 'start') % returns 8
%          idx = [2 3 8 9 10 11 33];
%          find_largemax(idx, '') % returns 10
%

function [idx_max, len] = find_largemax(idx, pos)
if ~isempty(idx)
    Mat(1, idx-min(idx)+1) = ones(1, length(idx));
    Diff = Mat(1:end-1) - Mat(2:end);
    Start = [0 find(Diff == -1)];
    Len = zeros(1, length(Start));
    Tmp = find(Diff(Start(1)+1:end) == 1, 1);
    if isempty(Tmp)
        Len(1) = length(Diff) - Start(1) + 1;
    else
        Len(1) = Tmp;
    end
    for i = 2:length(Start)
        Tmp = find(Diff(Start(i):end) == 1, 1);
        if isempty(Tmp)
            Len(i) = length(Diff) - Start(i) + 1;
        else
            Len(i) = Tmp - 1;
        end
    end
    len = max(Len);
    idx_max_len = find(Len == len);
    switch pos
        case 'start',
            idx_max = Start(idx_max_len(1)) + min(idx);
        case 'end',
            idx_max = Start(idx_max_len(1)) + Len(idx_max_len(1)) + min(idx) - 1;
        otherwise,
            idx_max = Start(idx_max_len(1)) + floor(Len(idx_max_len(1))/2) + min(idx);
    end
else
    idx_max = [];
    len = 0;
end
if isempty(idx_max)
    idx_max = NaN;
end
