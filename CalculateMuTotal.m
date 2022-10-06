function mu = CalculateMuTotal(formula,E)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
s = pwd;
cd('F:\matlab program\nist2matlab')
R = CalculateMassAC(formula,E);
mu  =R.TotalWithCoherent;
cd(s);
end

