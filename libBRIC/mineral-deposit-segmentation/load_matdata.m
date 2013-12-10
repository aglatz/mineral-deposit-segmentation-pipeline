function [Subjects, J, D, V, V_ref, Sens, Spec] = load_matdata(MatName, varargin)
load(MatName);
if exist('out', 'var')
    Ret = [out.Ret];
    Subjects = cell(length(Ret), 1);
    cnt = 1;
    for idx = 1:length(out)
        Tmp = out(idx).Subjects;
        Subjects(cnt:(cnt+length(Tmp)-1)) = Tmp;
        cnt = cnt + length(Tmp);
    end
end

if nargin > 1 && varargin{1}
    Ret = [Ret.CA];
end

J = [Ret.Jaccard];
D = [Ret.Dice];
V = [Ret.Vol];
V_ref = [Ret.Vol_ref];
Sens = [Ret.Sensitivity];
Spec = [Ret.Specificity];
