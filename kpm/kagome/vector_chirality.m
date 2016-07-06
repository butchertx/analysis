%Calculate chirality of a lattice of spins
%assume triangular lattice a1 = (1,0); a2 = (1/2, sqrt(3)/2)
tic
arg_list = argv ();
input_file = arg_list{1};
outfile = arg_list{2};
chifile = arg_list{3};
raw_spins = csvread(input_file);
spins = split_spins(raw_spins);%this is now a cell array of 3-vectors
N = length(spins)%number of sites
L = sqrt(N);%length of each dimension
a1 = [1, 0];
a2 = [-.5, sqrt(3)/2];
Sij = cell(L,L);
for j = 1:L
  for i = 1:L
    Sij{i,j} = spins{(i - 1)*L + j};
  end
end
chi = zeros(2*(L - 1),(L - 1));
for x = 1:(L - 1)
  for y = 1:(L - 1)
    %first triangle: x, x+1, y+1
    chi(2*x - 1, y) = dot(Sij{x, y+1}, cross(Sij{x, y}, Sij{x+1, y}));
    %second triangle: x+1, (x+1, y+1), y+1
    chi(2*x, y) = dot(Sij{x+1}, cross(Sij{x+1, y+1}, Sij{x, y+1}));
  end
end
plot_mat = zeros((L-1)*(L-1), 3);
for i = 1:(L-1)
  for j = 1:(L-1)
    plot_mat((i-1)*(L-1) + j, 1:2) = j*a1 + i*a2;
    plot_mat((i-1)*(L-1) + j, 3) = chi(2*i - 1,j);
  end
end
dims = size(chi);
avg_chi = sum(sum(chi))/(dims(1)*dims(2));
file = fopen(chifile);
fprintf(file,'%d\n',avg_chi);
save -ascii 'data.txt' plot_mat
toc
