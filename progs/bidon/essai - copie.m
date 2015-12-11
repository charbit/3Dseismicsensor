clear
% rotation axis eigenvector of R associated to the 
% eigenvalue 1
% angle de rotation theta verifies
%        Trace(R) = 1+2cos(theta)
%
au_deg =  [30;12;20];%rand(3,1)*180;%
au_rd = au_deg*pi/180;
cosau = cos(au_rd);
sinau = sin(au_rd);

eu_deg = rand(3,1)*180;%
eu_rd = eu_deg*pi/180;
coseu = cos(eu_rd);
sineu = sin(eu_rd);


Vu = [cosau(1)*coseu(1) sinau(1)*coseu(1) sineu(1); ...
    cosau(2)*coseu(2) sinau(2)*coseu(2) sineu(2); ...
    cosau(3)*coseu(3) sinau(3)*coseu(3) sineu(3)];


ar_deg =  [30;12;20];%rand(3,1)*180;%
ar_rd = ar_deg*pi/180;
cosar = cos(ar_rd);
sinar = sin(ar_rd);

er_deg = rand(3,1)*180;%
er_rd = er_deg*pi/180;
coser = cos(er_rd);
siner = sin(er_rd);


Vr = [cosar(1)*coser(1) sinar(1)*coser(1) siner(1); ...
    cosar(2)*coser(2) sinar(2)*coser(2) siner(2); ...
    cosar(3)*coser(3) sinar(3)*coser(3) siner(3)];

A = randn(3);
Suu = A*A';

