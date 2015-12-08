function V = matrixtrihadron(a_deg,e_deg)


a_rd  = a_deg*pi/180;
e_rd  = e_deg*pi/180;

cosa  = cos(a_rd);
sina  = sin(a_rd);

cose  = cos(e_rd);
sine  = sin(e_rd);


V    = zeros(3);
for k=1:3
V(k,:) = [cosa(k)*cose(k) sina(k)*cose(k) sine(k)];
end
