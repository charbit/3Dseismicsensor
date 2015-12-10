function voptim = extract1direction(...
    filteredsignalsUTk, filteredsignalsREF, ...
    HUTfk, HREFf, Vr)
%=======================================================================
% Synopsis:
%    voptim = extract1direction(...
%     filteredsignalsUTk, filteredsignalsREF, ...
%     HUTk, HREF, Vr)
% Inputs: 
%     - filteredsignalsUTk: 1 x N
%           signal of the selected SUT channel
%     - filteredsignalsREF: 3 x N
%           signal of the 3 SREF channels 
%     - HUTk: 1 x N
%           frequency response of the selected SUT channel
%     - HREF: 3 x N
%           frequency response of the 3 SREF channels
%     - Vr: unitay vectors of the 3 SREF channels
%           each row is an unitary vector.
% Output
%     - voptim: 3 x 1
%           unitary vector in the direction
%           of the selected SUT channel.
% 
%=======================================================================
XUTfk   = fft(filteredsignalsUTk,[],2);
XREFf   = fft(filteredsignalsREF,[],2);
Surk    = (XUTfk ./ HUTfk)  * XREFf';
Srr     = (XREFf ./ HREFf)  * XREFf';

A       = Srr'/(Vr');
ATA     = A'*A;
[UU,DD] = eig(ATA);

mu1     = DD(1,1);
mu2     = DD(2,2);
mu3     = DD(3,3);
W       = A'*(Surk');
gamma1  = abs(W'*UU(:,1))^2;
gamma2  = abs(W'*UU(:,2))^2;
gamma3  = abs(W'*UU(:,3))^2;

P(1)=1;
P(2)=2*mu1 + 2*mu2 + 2*mu3;
P(3)=mu1^2 + 4*mu1*mu2 + 4*mu1*mu3 + mu2^2 + ...
    4*mu2*mu3 + mu3^2 - gamma1 - gamma2 - gamma3;
P(4)=2*mu1*mu2^2 - 2*gamma2*mu1 - 2*gamma1*mu3 - ...
    2*gamma3*mu1 - 2*gamma2*mu3 - 2*gamma3*mu2 - ...
    2*gamma1*mu2 + 2*mu1^2*mu2 + 2*mu1*mu3^2 + ...
    2*mu1^2*mu3 + 2*mu2*mu3^2 + 2*mu2^2*mu3 + 8*mu1*mu2*mu3;
P(5)=mu1^2*mu2^2 - gamma2*mu1^2 - gamma1*mu3^2 - ...
    gamma3*mu1^2 - gamma2*mu3^2 - gamma3*mu2^2 - ...
    gamma1*mu2^2 + mu1^2*mu3^2 + mu2^2*mu3^2 + ...
    4*mu1*mu2*mu3^2 + 4*mu1*mu2^2*mu3 + 4*mu1^2*mu2*mu3 - ...
    4*gamma1*mu2*mu3 - 4*gamma2*mu1*mu3 - 4*gamma3*mu1*mu2;
P(6)=2*mu1^2*mu2^2*mu3 + 2*mu1^2*mu2*mu3^2 - ...
    2*gamma3*mu1^2*mu2 - 2*gamma2*mu1^2*mu3 + ...
    2*mu1*mu2^2*mu3^2 - 2*gamma3*mu1*mu2^2 - ...
    2*gamma2*mu1*mu3^2 - 2*gamma1*mu2^2*mu3 - ...
    2*gamma1*mu2*mu3^2;
P(7)=mu1^2*mu2^2*mu3^2 - gamma3*mu1^2*mu2^2 - ...
    gamma2*mu1^2*mu3^2 - gamma1*mu2^2*mu3^2;

valP=roots(P);
v=zeros(3,6);
J=zeros(6,1);
for ii=1:6
    v(:,ii) = inv(ATA+valP(ii)*eye(3))*W;
    J(ii) = norm(Surk-v(:,ii)'*A');
end

[bid,indmin]  = min(J);
voptim        = real(v(:,indmin));
%===========================================================