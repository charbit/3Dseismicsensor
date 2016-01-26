function vmin = extract1directionSCP(...
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

T = size(filteredsignalsUTk,2);
nbshift = 10;
wshift = T/nbshift;
Lw = wshift*2;
Suuf = 0;
MSC  = zeros(3,3,Lw);
Srrf  = zeros(3,3,Lw);
Sruf  = zeros(3,Lw);
winhan = hanning(Lw)';
winhan = winhan/ sqrt(sum(winhan .^2));
winhan3 = (ones(3,1)*winhan);
for iw=1:nbshift-1
    id1   = (iw-1)*wshift+1;
    id2   = id1 + Lw-1;
    xut   = filteredsignalsUTk(:,id1:id2) .* winhan;
    xrt   = filteredsignalsREF(:,id1:id2) .* winhan3;
    yuf   = fft(xut,[],2);% ./ HUTfk;
    yrf   = fft(xrt,[],2);% ./ HREFf;
    Suuf  = Suuf + abs(yuf) .^2;
    for ifq=1:Lw
        Srrf(:,:,ifq) = Srrf(:,:,ifq)+yrf(:,ifq)*yrf(:,ifq)';
        Sruf(:,ifq)   = Sruf(:,ifq)+yrf(:,ifq)*yuf(ifq)';
    end
end
% maxeig = zeros(Lw,1);
% for ifq=1:Lw
%     MSC(:,:,ifq)  = Srrf(:,:,ifq)\(Sruf(:,ifq)*Sruf(:,ifq)')/Suuf(:,ifq);
%     maxeig(ifq)   = max(eig(MSC(:,:,ifq)));
% end
% see document for notation
meanSrrf = mean(Srrf,3)/Lw;
meanSruf = mean(Sruf,2)/Lw;
Z       = meanSrrf'/(Vr');
ZHZ     = Z'*Z;
[UU,DD] = eig(ZHZ);
w       = Z'*(meanSruf);

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
    J_ii  = norm(meanSruf-Z*v_ii);
    if J_ii<Jmin
        Jmin = J_ii;
        vmin = real(v_ii);
    end
end

%===========================================================