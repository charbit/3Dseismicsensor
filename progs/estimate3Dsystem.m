clear
% attention the signal vectors and its Fourier transforms
% are ROW vetors.
rng(10);
%======= SREF
ar_deg    = [30;120;10];
er_deg    = [0.001; 0.002; 89.9];
Vr        = matrixtrihedron(ar_deg, er_deg);
%======= SUT randomly selected around the SREF
au_deg    = ar_deg+randn(3,1);
eu_deg    = er_deg+0.5*randn(3,1);
Vu        = matrixtrihedron(au_deg, eu_deg);
%==================================
% model of signals


% xu(t) = hu(t) * (Vu x g(t) + noise)
% xr(t) = hr(t) * (Vr x g(t) + noise)


%== generate N acceleration vectors
% and projections on SUT and SREF
N           = 1000;
wt          = randn(N,3);
[bnum,aden] = butter(4,2*[0.01 0.1]);
gt          = filter(bnum,aden,wt);
sut         = Vu*gt';
srt         = Vr*gt';
Pow         = mean([std(sut,[],2);std(sut,[],2)].^2);
%== we use SUT and SREF poles and zeroes, a littlebit random
% on poles, but you can also remove the randomness
sigmaPZ = 0.001;
zerosSUT      = [1;1;-1;-1]*ones(1,3);
polesSUT      = (1-sigmaPZ*rand)*...
    [0.9992+1j*0.0008;0.9992-1j*0.0008; -0.8841;-0.8699] * ...
    ones(1,3);
zerosSREF     = [1;1;-1;-1]*ones(1,3);
polesSREF     = (1-sigmaPZ*rand)*...
    [0.9992+1j*0.0008;0.9992-1j*0.0008; -0.8841;-0.8699] * ...
    ones(1,3);
%== all impulse responses vanish when time goes to infinity
% then we only consider a length of 100
hut = zeros(3,100);
hrt = zeros(3,100);
for ik=1:3
    numUT      = poly(zerosSUT(:,ik));
    denUT      = poly(polesSUT(:,ik));
    numREF     = poly(zerosSREF(:,ik));
    denREF     = poly(zerosSREF(:,ik));
    hut(ik,:)  = filter(numUT, denUT,eye(100,1));
    hrt(ik,:)  = filter(numREF, denREF,eye(100,1));
end
Huf    = fft(hut,N,2);
Hrf    = fft(hrt,N,2);
%===== Monte Carlo approach
Lruns  = 1000;
erra   = zeros(3,Lruns);
erre   = zeros(3,Lruns);
trihed = zeros(3,Lruns);
% noise level
SNRdB  = 20;
sigman = 10^(-SNRdB/20)*Pow;
for ir = 1:Lruns
    xutnoise = zeros(3,N);
    xrtnoise = zeros(3,N);
    for ik=1:3
        xutnoisebeforeF = sut(ik,:) + sigman*randn(1,N);
        xrtnoisebeforeF = srt(ik,:) + sigman*randn(1,N);
        xutnoise(ik,:)  = filter(hut(ik,:),1,xutnoisebeforeF);
        xrtnoise(ik,:)  = filter(hrt(ik,:),1,xrtnoisebeforeF);
    end
    %== start computation
    hatau_degk = zeros(3,1);
    hateu_degk = zeros(3,1);
    for kk=1:3
        voptim = extract1direction(...
            xutnoise(kk,:), xrtnoise, Huf(kk,:), Hrf, Vr);
        hatau_degk(kk) = real(atan2(voptim(2),voptim(1)))*180/pi;
        hateu_degk(kk) = real(asin(voptim(3)))*180/pi;
    end
    erra(:,ir)   = hatau_degk-au_deg;
    erre(:,ir)   = hateu_degk-eu_deg;
    trihed(:,ir) = voptim;
end

[mean(erra,2) std(erra,[],2)]
[mean(erre,2) std(erre,[],2)]

%%
figure(1)
for ii=1:3
    subplot(2,3,ii)
    DD = sprintf('D%i',ii);
    boxplot(erra(ii,:),'label',DD)
    set(gca,'fontname','times','fontsize',12)
    if ii==2
        ht = title(sprintf('the 3 azimut component\nerrors - degree'));
        pt = get(ht,'position');
        pt(2)=pt(2)*1.023;
        set(ht,'position',pt);
    end
    
    subplot(2,3,ii+3)
    boxplot(erre(ii,:),'label',DD);
    set(gca,'fontname','times','fontsize',12)
    if ii==2
        ht = title('the 3  component errors');
        ht = title(sprintf('the 3 elevation component\nerrors - degree'));
        pt = get(ht,'position');
        pt(2)=pt(2)*1.03;
        set(ht,'position',pt);
        textx = sprintf('%i samples, SNR = %i dB\nSimulation of %i runs',N,SNRdB,Lruns);
        hx = xlabel(textx);        
        pt = get(hx,'position');
        pt(2)=pt(2)*1.4;
        set(hx,'position',pt);

    end
end

HorizontalSize = 12;
VerticalSize   = 16;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');

printdirectory  = ' ../texte/';
fileprint = sprintf('%sboxplotoaande.eps',...
    printdirectory);

fileprintepscmd = sprintf('print -depsc -loose %s',fileprint);
fileeps2pdfcmd  = sprintf('!epstopdf %s',fileprint);
filermcmd       = sprintf('!rm %s',fileprint);


%     eval(fileprintepscmd)
%     eval(fileeps2pdfcmd)
%     eval(filermcmd)
