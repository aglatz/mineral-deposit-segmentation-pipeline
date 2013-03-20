function [RetAll, Over, Subjects] = segment_us_mp(SubjectFile, N_cpus)
% Function that parallelizes the segmentation (see pSegment.m) across
% CPU cores using pMatlab (see http://www.ll.mit.edu/mission/isr/pmatlab ).
% To circumvent the pMatlab installation execute pSegment_us.m directly.
% INPUTS: SubjectFile - A csv file in the format used for preprocessing. The
%                       last column is read and should contain the subject
%                       directory paths.
%         N_cpus - Number of CPU cores used for data processing.
% RETURNS: RetAll - A cell variable that contains limited mask statistics
%                   and validation results if a reference standard mask
%                   was found for the corresponding subject. See
%                   validate() for the format.
%          Over - Contains info about GRE cutoff thresholds calculated
%                 with the generated and reference masks, if given.
%                 See last two return values of segment_us().
%          Subjects - Cell array containing subject paths.

% Placeholder... will be filled in 'eval()' below.
RetAll = [];
Over = [];
Subjects = {};

% pMatlab
Cmd = sprintf('SubjectFile=''%s'';pSegment_us', SubjectFile);
eval(pRUN(Cmd, N_cpus, {}));
