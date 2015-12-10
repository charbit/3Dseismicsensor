function orientUT_deg = extractorient(...
    filteredsignalsUT3, filteredsignalsREF3, HUT3, HREF3,orientREF_deg, ...
    lista_deg,liste_deg)



ar_deg    = orientREF_deg.a;
er_deg    = orientREF_deg.e;
Vr        = matrixtrihadron(ar_deg, er_deg);
invVr     = inv(Vr);

nbvala    = length(lista_deg);
nbvale    = length(liste_deg);


XUTf  = fft(filteredsignalsUT3,[],2);
XREFf = fft(filteredsignalsREF3,[],2);
Sur   = (XUTf ./ HUT3)   * XREFf';
Srr   = (XREFf ./ HREF3) * XREFf';


J_ii  = zeros(nbvala,nbvale);
hatau_deg = zeros(3,1);
hateu_deg = zeros(3,1);

for ii=1:3
    Suu1D      = Sur(ii,:);
    for ia_deg = 1:nbvala
        a_deg  = lista_deg(ia_deg);
        for ie_deg = 1:nbvale
            e_deg  = liste_deg(ie_deg);
            V      = matrixtrihadron(a_deg, e_deg);
            J_ii(ia_deg,ie_deg) = norm(V*invVr*Srr-Suu1D);
        end
    end
    [ic,ir] = find(J_ii==min(min(J_ii)));
    hatau_deg(ii) = lista_deg(ic);
    hateu_deg(ii) = liste_deg(ir);
end
orientUT_deg.a = hatau_deg;
orientUT_deg.e = hateu_deg;

%===========================================================
function V = matrixtrihadron(a_deg,e_deg)

K     = length(a_deg);
a_rd  = a_deg*pi/180;
e_rd  = e_deg*pi/180;

cosa  = cos(a_rd);
sina  = sin(a_rd);

cose  = cos(e_rd);
sine  = sin(e_rd);


V     = zeros(K,3);
for k=1:K
    V(k,:) = [cosa(k)*cose(k) sina(k)*cose(k) sine(k)];
end
