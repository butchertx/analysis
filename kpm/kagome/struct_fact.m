function s_q = struct_fact(Lx, Ly, infile, plotfile)
%calculate the structure factor on the kagome lattice for the given file of
%input spins

b1 = 2*pi*[1, -1/sqrt(3)];%b1 dot a1 = 2*pi
b2 = 4*pi/sqrt(3)*[0, 1];%b2 dot a2 = 2*pi
basis1 = .5*[1, 0];
basis2 = .5*[.5, sqrt(3)/2];
%s(n,m) is at site n*a1 + m*a2, with three spins each

%indexed in kagome order, extract individual spin components from file
raw_spins = csvread(infile);
if(Lx*Ly*9 ~= length(raw_spins))
    Lx = sqrt(length(raw_spins)/9);
    Ly = Lx;
end
%make each site a cell array with components of each of 3 spins
spins = split_spins(raw_spins);

%now each lattice site has a matrix of spin components
for i = 1:length(spins)
    spins{i} = cell2mat(spins{i});
end
%spins is now a column of 1x9 vectors denoting 3 spins at each site
s_ij = cell(Lx, Ly);
for n = 1:Lx
    for m = 1:Ly
        s_ij{n, m} = spins{(m  - 1)* Lx + n};
    end
end

%put each component along 3rd dimension for fft
s_ij3D = zeros(Lx, Ly, 9);
for n = 1:Lx
    for m = 1:Ly
        s_ij3D(n, m, :) = s_ij{n, m};
    end
end

%transform a1 direction, then a2 direction
s_q3D = fft(fft(s_ij3D, [], 1), [], 2);

%exponential factors associated with basis vectors
expbasis1 = exp(1i*dot(basis1, b1));
expbasis2 = exp(1i*dot(basis2, b2));

%add the contributions from each sublattice
s_q_vec = s_q3D(:, :, 1:3) + s_q3D(:, :, 4:6)*expbasis1 + s_q3D(:, :, 7:9)*expbasis2;
s_q = zeros(Lx, Ly);


for n = 1:Lx
    for m = 1:Ly
        s_q(n, m) = real(dot(s_q_vec(n, m, :), s_q_vec(n, m, :)));
    end
end
s_q;

%s_q has each row denoting a row along b1. so X and Y will have the same
%number of rows and cols as s_q
interval1 = b1/Lx;
interval2 = b2/Ly;
figure()
hold on
colormap(redblue)

for n = 0:(Lx-1)
    for m = 0:(Ly-1)
        factor1 = m/Lx;
        factor2 = n/Ly;
        pos1 = factor1 * b1 + factor2 * b2;
        %4 cases: col outside BZ, row outside BZ, both outside, or none
        %outside
        %col outside BZ
%         if(factor2 >= .5*(1 + factor1) && factor2 >= 2*factor1)
%             pos1 = pos1 - b2;
%         %row outside BZ
%         elseif(factor2 <= .5*factor1 && factor2 <= 2*factor1 - 1)
%             pos1 = pos1 - b1;
%         %both outside
%         elseif(factor2 >= 1 - factor1 && factor2 >= .5*factor1 && factor2 <= 2*factor1)
%             pos1 = pos1 - b1 - b2;
%         end
        pos2 = pos1 + interval1;
        pos3 = pos2 + interval2;
        pos4 = pos3 - interval1;
        fill([pos1(1), pos2(1), pos3(1), pos4(1)],[pos1(2), pos2(2), pos3(2), pos4(2)], s_q(n + 1, m + 1))
    end
end

end