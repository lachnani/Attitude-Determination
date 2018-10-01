function [beta] = davenportq(vb,vn,w)
%DAVENPORTQ Determines attitude with Davenport's q method
%   Takes a cell array vb with body vectors, cell array vn with intertial
%   vectors, and w with weights for measurement accuracy. Outputs short
%   rotation quaternion with scalar as the fourth element. 

%Construct K matrix
n = size(vb);
b = w(1)*vb{1}*transpose(vn{1});
for i = 1:n
    b = b + (w(i)*vb{i}*transpose(vn{i}));
end
z = [b(2,3)-b(3,2);b(3,1)-b(1,3);b(1,2)-b(2,1)];
sig = trace(b);
s = b + transpose(b);
k = [sig z(1) z(2) z(3);z(1) s(1,1)-sig s(1,2) s(1,3);z(2) s(2,1) s(2,2)-sig s(2,3);z(3) s(3,1) s(3,2) s(3,3)-sig];
%Calculate and choose eigenvector
[vec,lam] = eig(k);
M = max(lam);
[M2,I] = max(M);
beta = vec(:,I);
%Make scalar fourth element
beta = circshift(beta,3);
%Make short rotation
if beta(4) < 0
    beta = -1*beta;
end
end


