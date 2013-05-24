function [RetAll,Subjects] = segment_us_mp(SubjectFile, RoiLabelTable, N_cpus, ...
                                           varargin)
% Function that calls pSegment.m (see its documentation for further info!).
% This function can parallelize the segmentation across CPU cores using
% pMatlab (see http://www.ll.mit.edu/mission/isr/pmatlab ) if N_cpus > 1.
% To circumvent the pMatlab installation execute pSegment_us.m directly.
% INPUTS: SubjectFile - A xls file in the format used for preprocessing. The
%                       last column is read and should contain the subject
%                       directory paths.
%         RoiLabelTable - A cell vector whose elements are label vectors.
%                         The first entry of the label vector are the
%                         labels of the ROIs that should be used to
%                         estimate the thresholds. All other ROI labels of
%                         the same vector should be segmented with the same
%                         threshold.
% OPTIONAL INPUTS: InterpFactor - A factor controlling the (fft)
%                                 interpolation. If =1 no
%                                 interpolation is done, >.1 entails
%                                 upsampling (to decrease PVE), <1 entails
%                                 downsampling.
%                  ReportName - File name where to store a report about the
%                               segmentation.
%                  ThreshFactor - A vector with two elements, the slope and
%                                 intercept, for changing the threshold.
%
%         N_cpus - Number of CPU cores used for data processing.
% RETURNS: RetAll - A cell variable that contains limited mask statistics
%                   and validation results if a reference standard mask
%                   was found for the corresponding subject. See
%                   validate() for the format.
%                 See last two return values of segment_us().
%          Subjects - Cell array containing subject paths.

% Placeholder... will be filled in 'eval()' below.
RetAll = [];
Over = [];
Subjects = {};

% Process optional args
N_vain = length(varargin);
for idx_vain = 1:N_vain
    arg_in = varargin{idx_vain};
    if ischar(arg_in) || isscalar(arg_in)
        switch arg_in
            case {'InterpFactor'}
                InterpFactor = varargin{idx_vain+1};
            case {'ReportName'}
                ReportName = varargin{idx_vain+1};
            case {'ThreshFactor'}
                ThreshFactor = varargin{idx_vain+1};
        end
    end
end

% Start processing
if N_cpus > 1
    % pMatlab
    InputFile = [tempname '.mat'];
    save(InputFile);
    idx = strfind(InputFile, '/');
    Cmd = sprintf('InputFile=''%s'';pSegment_us', InputFile(idx(end)+1:end));
    eval(pRUN(Cmd, N_cpus, {}));
    delete(InputFile);
else
    pSegment_us;
end
