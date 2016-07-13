%Run test functions from matlab window
defspins_infile = 'spins_in_test.csv';
defspins_outfile = 'spins_out_test.csv';
a1 = [1, 0];
a2 = [.5, sqrt(3)/2];
[Lx, Ly] = def_spins(defspins_infile, defspins_outfile);

raw_spins = csvread(defspins_outfile);
spins = split_spins(raw_spins);
Sij = cell(Lx,Ly);
for j = 1:Lx
  for i = 1:Ly
    Sij{i,j} = spins{(i - 1)*Lx + j};
  end
end

[x,y,z,u,v,w] = quiver_kagome_3D_planar(a1, a2, Sij, Lx, Ly);
quiver3(x,y,z,u,v,w,.5)

scalar_chirality(defspins_outfile, 'defspins_testplot', 'avg_chi.csv');

vector_chirality(defspins_outfile, 'defspins_vectestplot', 'avg_chi_vec.csv');