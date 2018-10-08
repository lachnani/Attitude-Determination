function [q] = olae(vb,vn,w)
%OLAE Determines attitude with OLAE method
%   Takes a cell array vb with body vectors, cell array vn with intertial
%   vectors, and w with weights for measurement accuracy. Outputs classical
%   Rodrigues parameter vector

%Constuct matrices
n = size(vb);
d = zeros(3*n(1),1);
s = zeros(3*n(1),3);
wm = zeros(3*n(1),3*n(1));
for i = 1:n(1)
    d(3*i-2:3*i) = vb{i}-vn{i};
    st = vb{i}+vn{i};
    s(3*i-2:3*i,:) = [0 -st(3) st(2) ; st(3) 0 -st(1) ; -st(2) st(1) 0 ];
    wm(3*i-2:3*i,3*i-2:3*i) = w(i)*eye(3);
end
%Calculate q
q = inv((transpose(s)*wm*s))*transpose(s)*wm*d;
end

