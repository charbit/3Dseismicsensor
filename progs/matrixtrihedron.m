function V = matrixtrihedron(a_deg,e_deg)
%=====================================================
% synopsis:
%     V = matrixtrihadron(a_deg,e_deg)
% V is a 3 x 3 matrix, whose rows are unitary vectors
% V(k,:) = [cosa(k)*cose(k) sina(k)*cose(k) sine(k)]
%=====================================================
%
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
%=====================================================
