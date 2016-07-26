function hamiltonian = band_structure(k)
%Test band structure code
N = 6;
S = 2;
J = .22*S;
a1 = [1, 0];
a2 = [.5, sqrt(3)/2];
a3 = a2 - a1;
Q1 = [2*pi, 0];
Q2 = [-pi, sqrt(3)*pi];
Q3 = [-pi, -sqrt(3)*pi];
hfind = @(k_)[0 cos(.5 * dot(k_,a2)) cos(.5 * dot(k_,a1))
    cos(.5 * dot(k_,a2)) 0 cos(.5 * dot(k_,a3))
    cos(.5 * dot(k_,a1)) cos(.5 * dot(k_,a3)) 0];
zeromat = zeros(N,N);
sigz = zeromat;
sigy = zeromat;
sigx = zeromat;

sigz(1,1) = 1;
sigz(4,4) = -1;
sigy(2,5) = -1i;
sigy(5,2) = 1i;
sigx(3,6) = 1;
sigx(6,3) = 1;
h_k0 = hfind(k);
h_kQ1 = hfind(k + Q1);
h_kQ2 = hfind(k + Q2);
h_kQ3 = hfind(k + Q3);
hmat = @(h_type) -2*[h_type zeros(3,3)
    zeros(3,3) h_type];

hamiltonian = [hmat(h_kQ3) zeromat zeromat -J*sigz
    zeromat hmat(h_kQ2) zeromat -J*sigy
    zeromat zeromat hmat(h_kQ1) -J*sigx
    -J*sigz -J*sigy -J*sigx hmat(h_k0)];
    
end