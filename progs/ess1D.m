clear

%======= gridding the 2 parameters
nbvala      = 300;
nbvale      = 320;
lista_deg   = linspace(0,360,nbvala);
liste_deg   = linspace(-90,90,nbvale);

%======= SUT
au_deg =  [230;81;50];%[360*rand(2,1);60];
eu_deg =  [43; 60;35];%[180*rand(2,1);35];
Vu     = matrixtrihadron(au_deg, eu_deg);

%======= SREF random
ar_deg    = 360*rand(3,1);%
er_deg    = 180*rand(3,1);%
Vr        = matrixtrihadron(ar_deg, er_deg);
invVr     = inv(Vr);
orientREF_deg.a = ar_deg  ;
orientREF_deg.e = er_deg  ;
%==================================
% model of signals
% Xu(f) = Hu(f) * Vu * G(f)
% Xr(f) = Hr(f) * Vr * G(f)
% Thanks to the responses Hu et Hr 
% which are assumed te be known
% we simulate Xuf and Xrf using random frequency responses.
% 

N           = 1000;
wt          = randn(N,3);
[bnum,aden] = butter(4,2*[0.01 0.1]);
gt          = filter(bnum,aden,wt);
gt          = gt';
Gf          = fft(gt,[],2);

hu          = randn(3,N)/N;
Huf         = fft(hu,[],2);
hr          = randn(3,N)/N;
Hrf         = fft(hr,[],2);

sigman      = 0.002;
Xuf         = (Huf .* (Vu * Gf)) + sigman*(randn(3,N)+1j*randn(3,N));
Xref        = (Hrf .* (Vr * Gf)) + sigman*(randn(3,N)+1j*randn(3,N));

xut         = real(ifft(Xuf,[],2));
xref        = real(ifft(Xref,[],2));

orientUT_deg = extractorient(xut,xref, ...
    Huf, Hrf ,orientREF_deg, ...
    lista_deg,liste_deg)
hatau_deg = orientUT_deg.a;
hateu_deg = orientUT_deg.e;
[hatau_deg au_deg]
[hateu_deg eu_deg]

