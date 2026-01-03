function F = root_b(b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global k_ini n ratio
F = -((b^(1 + 2*k_ini) + k_ini*b^k_ini - k_ini*b^(2 + k_ini) + b^(1 + n)*(-1 + n)-...
  n*b^n)/((-1 + b)*(b^k_ini - b^(2*k_ini) + b^(1 + k_ini) - b^n)))-ratio;
end

