function vmin = extract1direction(...
    filteredsignalsUTk, filteredsignalsREF, ...
    HUTfk, HREFf, Vr)
%=======================================================================
% Synopsis:
%    vmin = extract1direction(...
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
%     - vmin: 3 x 1
%           unitary vector in the direction
%           of the selected SUT channel solution 
%           of the least square problem.
% Used functions:
%      - roots, poly, norm from Matlab
%
%=======================================================================
XUTfk   = fft(filteredsignalsUTk,[],2);
XREFf   = fft(filteredsignalsREF,[],2);
yurk    = (XUTfk ./ HUTfk) ; % deconvolution
yREFf   = (XREFf ./ HREFf) ; % deconvolution
% see document for notation
Z       = yREFf'/(Vr');
ZHZ     = Z'*Z;
[UU,DD] = eig(ZHZ);
w       = Z'*(yurk');

mu1     = DD(1,1);
mu2     = DD(2,2);
mu3     = DD(3,3);
gamma1  = abs(w'*UU(:,1))^2;
gamma2  = abs(w'*UU(:,2))^2;
gamma3  = abs(w'*UU(:,3))^2;

Polylambda    = zeros(7,1);
Polylambda(1) = 1;
Polylambda(2) = 2*mu1 + 2*mu2 + 2*mu3;
Polylambda(3) = mu1^2 + 4*mu1*mu2 + 4*mu1*mu3 + mu2^2 + ...
    4*mu2*mu3 + mu3^2 - gamma1 - gamma2 - gamma3;
Polylambda(4) = 2*mu1*mu2^2 - 2*gamma2*mu1 - 2*gamma1*mu3 - ...
    2*gamma3*mu1 - 2*gamma2*mu3 - 2*gamma3*mu2 - ...
    2*gamma1*mu2 + 2*mu1^2*mu2 + 2*mu1*mu3^2 + ...
    2*mu1^2*mu3 + 2*mu2*mu3^2 + 2*mu2^2*mu3 + 8*mu1*mu2*mu3;
Polylambda(5)=mu1^2*mu2^2 - gamma2*mu1^2 - gamma1*mu3^2 - ...
    gamma3*mu1^2 - gamma2*mu3^2 - gamma3*mu2^2 - ...
    gamma1*mu2^2 + mu1^2*mu3^2 + mu2^2*mu3^2 + ...
    4*mu1*mu2*mu3^2 + 4*mu1*mu2^2*mu3 + 4*mu1^2*mu2*mu3 - ...
    4*gamma1*mu2*mu3 - 4*gamma2*mu1*mu3 - 4*gamma3*mu1*mu2;
Polylambda(6)=2*mu1^2*mu2^2*mu3 + 2*mu1^2*mu2*mu3^2 - ...
    2*gamma3*mu1^2*mu2 - 2*gamma2*mu1^2*mu3 + ...
    2*mu1*mu2^2*mu3^2 - 2*gamma3*mu1*mu2^2 - ...
    2*gamma2*mu1*mu3^2 - 2*gamma1*mu2^2*mu3 - ...
    2*gamma1*mu2*mu3^2;
Polylambda(7)=mu1^2*mu2^2*mu3^2 - gamma3*mu1^2*mu2^2 - ...
    gamma2*mu1^2*mu3^2 - gamma1*mu2^2*mu3^2;

nbroots = 6;
rootsP  = roots(Polylambda);
Jmin    = inf;
vmin    = NaN;
for ii=1:nbroots
    v_ii  = (ZHZ+rootsP(ii)*eye(3))\w;
    J_ii  = norm(yurk-v_ii'*Z');
    if J_ii<Jmin
        Jmin = J_ii;
        vmin = real(v_ii);
    end
end

%===========================================================