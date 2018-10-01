function [q] = quest(vb,vn,w,iter)
%QUEST Determines attitude with QUEST method
%   Takes a cell array vb with body vectors, cell array vn with intertial
%   vectors, and w with weights for measurement accuracy. Outputs classical
%   Rodrigues parameter vector

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
%Estimate maximum eigenvalue
lam = sum(w);
f = det(k-lam*eye(4));
fp = trace(-1*f./(k-lam*eye(4)));
for j = 1:iter
    lam = lam - (f/fp);
    f = det(k-lam*eye(4));
    fp = trace(-1*f./(k-lam*eye(4)));
end
%Calculate CRPs
q = z/((lam+sig)*eye(3));
end
%Code not finished
