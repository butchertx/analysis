function [Lx, Ly] = def_spins(infile, outfile)
%outfile is the desired name for the output file (a <spins>.csv-type file)
%infile is the input file defining the desired spin structure
inputmat = csvread(infile, 0, 1);

dim = inputmat(2, :);
Lx = dim(1);
Ly = dim(2);
delta_sub1 = inputmat(4:6, :);
delta_sub2 = inputmat(8:10, :);
delta_sub3 = inputmat(12:14, :);
spins = cell(1, Lx * Ly);
zetas = zeros(1, 3); %cos(Q_eta dot r)

for l = 0:(Lx - 1)
    zetas(3) = 1 - 2 * mod(l, 2);
    for n = 0:(Ly - 1)
        zetas(1) = 1 - 2 * mod(n, 2);
        if(zetas(1) == zetas(3))
            zetas(2) = 1;
        else
            zetas(2) = -1;
        end
        spin1 = delta_sub1(1,:) * zetas(1) + delta_sub1(2,:) * zetas(2) + delta_sub1(3,:) * zetas(3);
        spin2 = delta_sub2(1,:) * zetas(1) + delta_sub2(2,:) * zetas(2) + delta_sub2(3,:) * zetas(3);
        spin3 = delta_sub3(1,:) * zetas(1) + delta_sub3(2,:) * zetas(2) + delta_sub3(3,:) * zetas(3);
        spins{n * Lx + l + 1} = [spin1, spin2, spin3];
    end
end
spins = cell2mat(spins);
csvwrite(outfile,spins);
end