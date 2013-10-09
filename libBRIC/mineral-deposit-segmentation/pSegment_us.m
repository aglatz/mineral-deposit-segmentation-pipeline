% This script gets called from segment_us_mp() and distributes the
% input data processing with segment_us_single() or function with
% the same input/output arguments accross several CPU cores using
% pMatlab (see http://www.ll.mit.edu/mission/isr/pmatlab/ for further info).
%

% First include our libraries
addpath('../misc-matlab/');

% Include external libraries - create symbolic links if it fails here!
addpath('NIFTI/');
addpath('LIBRA/');

close all; % No 'clear all' otherwise we loose our input variables!

% Load input vars in mp mode
InputFileVar = 'InputFile';
if exist(InputFileVar, 'var') > 0
    PMatDir = tempdir;
    InputFile = fullfile(PMatDir, InputFile);
    load(InputFile);
end

% Default input arguments
if exist('ReportName', 'var') <= 0 || isempty(ReportName)
    ReportName = [];
end
if exist('InterpFactor', 'var') <= 0 || isempty(InterpFactor)
    InterpFactor = 1;
end
if exist('ThreshFactor', 'var') <= 0 || isempty(ThreshFactor)
    ThreshFactor = [1 0];
end
if exist('AdaptiveFlag', 'var') <= 0 || isempty(AdaptiveFlag)
    AdaptiveFlag = true;
end
if exist('IntvarP', 'var') <= 0 || isempty(IntvarP)
    IntvarP = 0.5;
end
if exist('FuncName', 'var') <= 0 || isempty(FuncName)
    FuncName = 'segment_us_single';
end

% Read subject file - last column contains the paths to the subject dirs
[Subjects, N_total] = get_subject_dirs(SubjectFile);

if exist(InputFileVar, 'var') > 0
    % pMatlab setup - Init shared matrix
    Wmap=map([Np 1], {}, 0:Np-1);
    % Wmap=1; % Uncomment for sequential processing
    Over = zeros(N_total, length(RoiLabelTable)*2, Wmap);
    % Get portion of matrix for this process
    OverLoc = local(Over);
    IdxLoc = global_ind(Over, 1);
else
    OverLoc = zeros(N_total, length(RoiLabelTable)*2, 1);
    IdxLoc = 1:N_total;
end
for idx_j = 1:size(OverLoc, 1);
	idx_sub = IdxLoc(idx_j);
	Subject = Subjects{idx_sub};
    fprintf('Subject: %s ...\n', Subject);
    
	fh = str2func(FuncName);
	Ret = fh(Subject, RoiLabelTable, ReportName, InterpFactor, ...
             ThreshFactor, AdaptiveFlag, IntvarP);
	save(fullfile(Subject, 'Ret.mat'), 'Ret');
end

if exist(InputFileVar, 'var') > 0
    % pMatlab exit - Aggregate shared matrix
    Over = put_local(Over, OverLoc);
    Over = agg(Over);
else
    Over = OverLoc;
end

% Collect result data in parent process (Pid=0)
if ~exist('Pid','var') || ~Pid
    % Init return structure
    RetAll = Ret(end);
    % Combine all returned structures to a array of structures
    for i = 1:N_total
        RetPath = fullfile(char(Subjects(i)), 'Ret.mat');
        load(RetPath);
        RetAll(i) = Ret(end); % Last one is best one
    end
end
