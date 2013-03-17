function [V] = get_volume(SM, varargin)
% Returns volume of given mask.
% INPUTS: SM - Binary mask
%         varargin{1} - voxel size
% OUTPUT: V - volume
%
if nargin > 1
    F = varargin{1};
else
    F = [1 1 1];
end
V = sum(sum(sum(SM))) * prod(F);

