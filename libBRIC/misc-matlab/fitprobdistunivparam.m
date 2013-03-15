function [Bfit, Bfit_hist] = fitprobdistunivparam(X, Dist_type, N_min, Tf)
% Scales the data in input cell vector 'X' to conver the range 'Tf' unless the
% vector contains less than 'N_min' elements. Then uses Matlab's fitdist()
% function to fit specified distributions (see help fitdist) to the data.
% INPUTS: X - Cell vector, of which each element contains an input data vector
%         Dist_type - Cell array that contains valid distribution
%                     types as string (see help fitdist)
%         N_min - Minimum length of each data vector
%         Tf - Scaling range (e.g. [0 1] for fitting a beta distribution)
% OUTPUTS: Bfit - Cell matrix that has the same length as 'X' and of which
%                 each row contains the best fitting distribution (according
%                 to the Akakine Information Criteria)
%          Bfit_hist - Histogram which shows the frequency with which each
%                      distribution of 'Dist_type' was chosen as best
%                      fitting distribution.
%
N = length(X);
Bfit = cell(N, 3);
N_dist = length(Dist_type);
for idx=1:N
    if length(X{idx}) >= N_min % Everything above minimum sample size
                               % is OK for fitting
        % Linearily transform data so it lies in the interval Tf
        X_dbl = double(X{idx});
        P = polyfit([min(X_dbl) max(X_dbl)], Tf, 1);
        X_tf = polyval(P, X_dbl);
        
        % Now fit the distribution to each data element of the
        % input cell vector
        PD = cell(N_dist, 1);
        Res = zeros(N_dist, 1);
        for idx_dist = 1:N_dist
            PD{idx_dist} = fitdist(X_tf(:), Dist_type{idx_dist});
            k = numel(PD{idx_dist}.Params);
            NLL = PD{idx_dist}.NLogL;
            Res(idx_dist) = -2*(-NLL)+2*k; % AIC: Akakine Information Criteria
                                           % for comparing the models
        end
        [Tmp, Idx] = sort(Res);
        Bfit{idx, 1} = Idx(1); % Model with the lowest AIC is best
        Bfit{idx, 2} = PD{Bfit{idx, 1}}; % Save distribution object
        Bfit{idx, 3} = length(X_dbl); % store the number of samples
    else
        Bfit{idx, 1} = -1;
    end
end
% Calculate histogram
Bfit_hist = hist(cell2mat(Bfit(:, 1)), [-1 1:N_dist]) ./ N .* 100;
