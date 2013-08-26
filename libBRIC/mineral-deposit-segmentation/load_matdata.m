function [Subjects, J, D, V, V_ref, IntVarP] = load_matdata(MatName)
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

J = [Ret.Jaccard];
D = [Ret.Dice];
V = [Ret.Vol];
V_ref = [Ret.Vol_ref];
Input = [Ret.Input];
IntVarP = [Input.IntvarP];
