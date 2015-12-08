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
