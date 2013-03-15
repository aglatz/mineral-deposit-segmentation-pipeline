function [Lab, Loc] = conncomp_mask(L, S_roi, varargin)
% Determines which ROI in 'S_roi' mostly contains the
% connected components in 'L'.
% INPUTS: L - the connected component mask as returned
%             by conncomp_init().
%         S_roi - the ROI mask where each ROI has a
%                 unique label.
% OPTIONAL INPUTS: varargin{1} - the fraction that a
%                  connected component has to intersect
%                  with a ROI to be counted part of it.
%                  varargin{2} - A volume with weights
%                  that is used to resolve conflicts
%                  when a connected component intersects
%                  with two or more ROIs equally.
% RETURNS: Lab - List of connected component labels
%          Loc - List of locations (ROI labels)
%
outside = 1; % Any region which intersects with ROI
if nargin > 2
    if varargin{1} < 0 || varargin{1} > 1
        fprintf('Range error for outside! Setting to 1.\n');
    else
        outside = varargin{1};
    end
end
if nargin > 3
    S_weight = varargin{2};
else
    S_weight = ones(size(L));
end

Lab = unique(L)';
Lab = Lab(2:end); % Exclude background
N_lab = length(Lab);
Lab_roi = unique(S_roi)';
Lab_roi = Lab_roi(2:end); % Exclude background
N_lab_roi = length(Lab_roi);
SM_roi = logical(S_roi);
Loc = ones(size(Lab)) .* -1; % -1 means outside ROI
for idx = 1:N_lab
    SM = L == Lab(idx);    
    % Is inside ROI?
    V_in = get_volume(SM_roi & SM);
	V_out = get_volume(~SM_roi & SM);
    if V_out/(V_out+V_in) < outside
        % Which region does SM belong to?
        S = single(S_roi .* cast(SM, class(S_roi)));
        Vol = hist(S(:), single([0 Lab_roi]));
        idx_max = find(Vol(2:end) == max(Vol(2:end)));
        if length(idx_max) > 1
            % Use weight matrix to resolve
            I = nan(size(Lab_roi));
            for idx_roi = 1:N_lab_roi
                SM_tmp = S == Lab_roi(idx_roi);
                if sum(SM_tmp(:)) > 0
                    I(idx_roi) = median(single(S_weight(SM_tmp)));
                end
            end
            % We are searching the most hypointense part
            idx_I_max = find(I(idx_max) == min(I(idx_max)));
            if length(idx_I_max) == 1
                Loc(idx) = Lab_roi(idx_max(idx_I_max));
            else
                % Just pick first best
                Loc(idx) = Lab_roi(idx_max(1));
            end
        else
            % Simple case
            Loc(idx) = Lab_roi(idx_max);
        end
    end
end


