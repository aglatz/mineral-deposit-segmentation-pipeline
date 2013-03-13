function [OccCorr, BinCentre, X, Y] = plot_hist(I, I_Max, I_Max_Pos, ...
												hist_colour, print_kde)
% Generates and plots a histogram with optimal bin width and adds
% the kernel density estimate, if requested.
% INPUTS: I - the data vector to plot
%         I_Max - Used for scaling the y-axis of the histogram if not empty
%         I_Max_Pos - X-axis position of the maximum
%         hist_colour - Colour of the histogram
%         print_kde - flag that indicates whether a kernel density estimate
%                     should be added to the histogram
% RETURNS: OccCorr - Corrected Occurrence (y values of histogram bars)
%          BinCentre - Center values of the histogram bars
%          X - x values where the kernel density estimates were caclculated
%          Y - kernel density estimates ordered according to corresponding
%              x values.
%
I = double(I(:)); % hist() either wants float or double and a vector
len = length(I);
if len > 1
    len_thr = 30; % Empiric...
    nBins = [];
    if len > len_thr
        nBins = sshist(I);
    end
    if isempty(nBins)
        nBins = 15;
    end
    [Occ, BinCentre] = hist(I, nBins);

    if ~isempty(I_Max)
        MaxOcc = interp1(BinCentre, Occ, I_Max_Pos);
        OccCorr = Occ / MaxOcc * I_Max; % Scale to given max value
    else
        OccCorr = Occ;
    end

    if ~isempty(hist_colour)
        bar(BinCentre, OccCorr, hist_colour); % Histogram
    else
        bar(BinCentre, OccCorr); % Histogram
    end

    I(isnan(I)) = []; % kde() doesn't support NaNs
    if len > len_thr && sum(I) ~= 0 && print_kde
        hold on;
        [tmp1, Y_tmp, X] = kde(I, 2^14); clear tmp1;
        Y = Y_tmp / max(Y_tmp) * max(OccCorr); % Scale to given max value
        plot(X, Y, 'g');
    else
        X = [];
        Y = [];
    end
else
    OccCorr = [];
    BinCentre = [];
    X = [];
    Y = [];
end

