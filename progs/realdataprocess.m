clear
% attention the signal vectors and its Fourier transforms
% are ROW vetors.

load('year2015month10day22.mat')

[num,den] = butter(2,2*[0.01 1/4]);
%======= SREF

% theoretically if azimutREF_deg-elevationREF_deg 
% stays equal to +/-90°, results are unchanged
Ox.azimutREF_deg    = 0;
Ox.elevationREF_deg = 0;

Oy.azimutREF_deg    = 90;
Oy.elevationREF_deg = 0;

% theoretically the Oz azimut has no effect 
% because the Oz elevation is 90°
Oz.azimutREF_deg    = 0;
Oz.elevationREF_deg = 90;

azimutREF_deg       = [Ox.azimutREF_deg;...
    Oy.azimutREF_deg;Oz.azimutREF_deg];
elevationREF_deg    = [Ox.elevationREF_deg; ...
    Oy.elevationREF_deg; Oz.elevationREF_deg];
VREF = matrixtrihedron(azimutREF_deg, elevationREF_deg);
grammianVREF = VREF*VREF';
%==================================
N_sec     = 1000;
T         = size(signals_centered,1);
N         = fix(Fs_Hz*N_sec);
Nshift    = fix(N*0.5);
Lruns     = fix(T/Nshift);
Huf       = ones(3,N);
Hrf       = ones(3,N);

signals_filtered = zeros(size(signals_centered));
for ik=1:6
    signals_filtered(:,ik) = filter(num,den,signals_centered(:,ik));
end

hatazimutSUT_deg       = zeros(3,Lruns);
hatelevationSUT_deg    = zeros(3,Lruns);
for ir=1:Lruns
    id1                = (ir-1)*Nshift+1;
    id2                = min([id1+N-1,T]);   
    N_ir               = id2-id1+1;
    xutnoise           = zeros(3,N_ir);
    xrtnoise           = zeros(3,N_ir);
    xutnoise(1,1:N_ir) = signals_filtered(id1:id2,1);
    xutnoise(2,1:N_ir) = signals_filtered(id1:id2,2);
    xutnoise(3,1:N_ir) = signals_filtered(id1:id2,3);
    xrtnoise(1,1:N_ir) = signals_filtered(id1:id2,4);
    xrtnoise(2,1:N_ir) = signals_filtered(id1:id2,5);
    xrtnoise(3,1:N_ir) = signals_filtered(id1:id2,6);
    
    %== start computation
    for kk=1:3
        voptim = extract1direction(...
            xutnoise(kk,:), xrtnoise, Huf(kk,1:N_ir), Hrf(:,1:N_ir), VREF);
        hatazimutSUT_deg(kk,ir) = real(atan2(voptim(2),voptim(1)))*180/pi;
        hatelevationSUT_deg(kk,ir) = real(asin(voptim(3)))*180/pi;
    end
end
%%
figure(1)
for kk=1:3
    subplot(2,3,kk)
    hist(hatazimutSUT_deg(kk,:),10)
    subplot(2,3,kk+3)
    hist(hatelevationSUT_deg(kk,:),10)
end

meanazSUT = mean(hatazimutSUT_deg,2);
stdazSUT  = std(hatazimutSUT_deg,[],2);

meanelSUT = mean(hatelevationSUT_deg,2);
stdelSUT  = std(hatelevationSUT_deg,[],2);

hatVSUT = matrixtrihedron(meanazSUT,meanelSUT);

grammianVSUT = hatVSUT*hatVSUT';

grammianVSUT
grammianVREF
