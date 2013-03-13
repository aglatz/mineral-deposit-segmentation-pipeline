function [Ret] = get_appearance(I_roi, I_fe, P_fe, a)
% Returns a number (0,1,2) that corresponds to the
% appearance of the lesion w/r to the normal 
% tissue intensity. A ranksum() test is performed first.
% (see combine_subject_data script).
% INPUTS: I_roi - mean intensity of the normal-appearing tissue
%         I_fe - mean intensity of the lesion
%         P_fe - p-value from the ranksum test
%         a - significance level
% RETURNS: Ret - 0,1,2 acording to appearance
%
Ret = NaN; % N/A
if ~isnan(I_fe) && ~isnan(I_roi)
    if P_fe < a
        if I_fe < I_roi
            Ret = 1; % hypointense
        else
            Ret = 2; % hyperintense
        end
    else
        Ret = 0; % isointense
    end
end
