clear
% rotation axis eigenvector of R associated to the 
% eigenvalue 1
% angle de rotation theta verifies
%        Trace(R) = 1+2cos(theta)
%
a_deg =  [30;12;20];%rand(3,1)*180;%
a_rd = a_deg*pi/180;
cosa = cos(a_rd);
sina = sin(a_rd);
R1 = [1 0 0;0 cosa(1) -sina(1);0 sina(1) cosa(1)];
R2 = [cosa(2) 0 sina(2);0 1 0; -sina(2) 0 cosa(2)];
R3 = [cosa(3) -sina(3) 0;sina(3) cosa(3) 0;0 0 1];
R  = R1*R2*R3;
N  = 1000;
eig(R)
% Xu(f) = R * Hu(f) * s(f)
% Xr(f) = Hr(f) * s(f)
% Thanks to the responses Hu et Hr assumed te be known
% we have
%      s(f) =  (Hr^T(f)Hr(f))^-1Hr^T(f)Xr(f)
%
% ==>   Xu(f) = R * Hu(f) * (Hr^T(f)Hr(f))^-1Hr^T(f)Xr(f)

sf = randn(N,3)+1j*randn(N,3);
zf = sf * R;

