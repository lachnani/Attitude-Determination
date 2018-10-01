function [bn] = triad(b1,b2,i1,i2)
%TRIAD Estimates attitude with the TRIAD method
%   Takes, in order, 1st vector in body frame, 2nd vector in body frame,
%   1st vector in inertial frame, and 2nd vector in intertial frame.
%   Outputs inertial to body DCM. All inputs must be column vectors.

%Normalize vectors
b1 = b1/norm(b1);
b2 = b2/norm(b2);
n1 = i1/norm(i1);
n2 = i2/norm(i2);
%BT matrix
bt2 = cross(b1,b2)/norm(cross(b1,b2));
bt3 = cross(b1,bt2);
bt = [b1 bt2 bt3];
%NT matrix
nt2 = cross(n1,n2)/norm(cross(n1,n2));
nt3 = cross(n1,nt2);
nt = [n1 nt2 nt3];
%Calculate BN
bn = bt*transpose(nt);
end

