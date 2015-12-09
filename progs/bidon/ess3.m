clear
%======= SUT
au_deg =  180*rand(3,1);
eu_deg =  180*rand(3,1);
Vu     = matrixtrihadron(au_deg, eu_deg);

%======= SREF
ar_deg =  180*rand(3,1);
er_deg =  180*rand(3,1);
Vr     = matrixtrihadron(ar_deg, er_deg);


%==================================================================
% model of signals
% Xu(f) = Hu(f) * Vu * G(f)
% Xr(f) = Hr(f) * Vr * G(f)
% Thanks to the responses Hu et Hr assumed te be known
% we can assume Identity

N = 100;
Gf = randn(3,N)+1j*randn(3,N);
Xuf = Vu * Gf;
Xrf = Vr * Gf;

nbval = 10;
lista_deg = linspace(0,180,nbval);
liste_deg = linspace(0,180,nbval);
J = zeros(nbval,nbval,nbval,nbval,nbval,nbval);
for ia_deg1 = 1:nbval
    a_deg1 = lista_deg(ia_deg1);
    for ia_deg2 = 1:nbval
        a_deg2 = lista_deg(ia_deg2);
        for ia_deg3 = 1:nbval
            a_deg3 = lista_deg(ia_deg3);
            a_deg  = [a_deg1 a_deg2 a_deg3];
            for ie_deg1 = 1:nbval
                e_deg1 = liste_deg(ie_deg1);
                for ie_deg2 = 1:nbval
                    e_deg2 = liste_deg(ie_deg2);
                    for ie_deg3 = 1:nbval
                        e_deg3 = liste_deg(ie_deg3);
                        e_deg  = [e_deg1 e_deg2 e_deg3];
                        V = matrixtrihadron(a_deg, e_deg);
                        J(ia_deg1,ia_deg2,ia_deg3,ie_deg1,ie_deg2,ie_deg3) = norm(V*inv(Vr)*Xrf-Xuf);
                    end
                end
            end
        end
    end
end

                %
                %
                % zf = sf * R;
                %
