function [SubjectDirs, N_total] = get_subject_dirs(SubjectFile)
% Returns a cell array with the paths to the individual subject dirs.
% INPUTS: SubjectFile - Either xls or csv file with subject dirs
%                       in last column
% OUTPUTS: SubjectDirs - Cell array with subject dirs
%          N_total - number of subjects
%
N_total = 0;
try
    [~, ~, raw] = xlsread(SubjectFile); clear ndata text;
    
    % Get subject directories from last column
    SubjectDirs = raw(:, size(raw, 2));
    N_total = length(SubjectDirs);
catch
	fd = fopen(SubjectFile);
    if fd >= 0
        % Read all lines, ignore comments
        raw = textscan(fd,	'%s', ...
                            'delimiter', '\r\n', ...
                            'commentStyle', '#');
        fclose(fd);
        raw = raw{1};
        N_lines = length(raw);
        
        % Get subject directories from last column
        SubjectDirs = cell(N_lines, 1);
        is_first = true;
        N_total = 0;
        for idx_line = 1:N_lines
            line = textscan(raw{idx_line}, '%q', ...
                                           'delimiter', ',');
            line = line{1};
            if is_first
                N_elem = length(line);
                is_first = false;
            else
                if length(line) ~= N_elem
                    continue; % skip odd lines
                end
            end
            SubjectDirs(idx_line) = line(end);
            N_total = N_total + 1;
        end
    end 
end
if N_total <= 0
    error('pSegment_us:open', 'Cannot open input file!');
end
