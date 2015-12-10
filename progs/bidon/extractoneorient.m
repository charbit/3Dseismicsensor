function orientUTk_deg = extractoneorient(...
    filteredsignalsUTk, filteredsignalsREF3, HUTk, ...
    HREF3D,orientREF_deg, ...
    lista_deg,liste_deg)



ar_deg    = orientREF_deg.a;
er_deg    = orientREF_deg.e;
Vr        = matrixtrihadron(ar_deg, er_deg);
invVr     = inv(Vr);

nbvala    = length(lista_deg);
nbvale    = length(liste_deg);

XUTkf     = fft(filteredsignalsUTk,[],2);
XREFf     = fft(filteredsignalsREF3,[],2);
Surk      = (XUTkf ./ HUTk)   * XREFf';
Srr       = (XREFf ./ HREF3D) * XREFf';

Jk        = zeros(nbvala,nbvale);

for ia_deg = 1:nbvala
    a_deg  = lista_deg(ia_deg);
    for ie_deg = 1:nbvale
        e_deg  = liste_deg(ie_deg);
        V      = matrixtrihadron(a_deg, e_deg);
        Jk(ia_deg,ie_deg) = norm(V*invVr*Srr-Surk);
    end
end
[ic,ir] = find(Jk==min(min(Jk)));
orientUTk_deg.a = lista_deg(ic);
orientUTk_deg.e = liste_deg(ir);

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
