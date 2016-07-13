function [X,Y,Z,U,V,W] = quiver_kagome_3D_planar(a1, a2, Sij, Lx, Ly)
%take the two lattice vectors a1,a2 and the cell array of spins Sij
%at each point i*a2 + j*a1, Sij contains 3 3-vectors -- only need 
%the xy projection (first two components) and will assume basis
%vectors at each site to be .5*a1 and .5*a2

X = cell(size(Sij));
Y = cell(size(Sij));
Z = cell(size(Sij));
U = cell(size(Sij));
V = cell(size(Sij));
W = cell(size(Sij));

for i = 1:Ly
	for j = 1:Lx
		triangle = Sij{i,j};%cell array with 3 3-vectors
		position0 = (i-1)*a2 + (j-1)*a1;%n=2 vector
		position1 = position0 + .5*a1;
		position2 = position0 + .5*a2;
		
		X{i, j} = [position0(1), position1(1), position2(1)];
		Y{i, j} = [position0(2), position1(2), position2(2)];
        Z{i, j} = [0, 0, 0];

		triangle0 = triangle{1};%n=3 vector
		triangle1 = triangle{2};
		triangle2 = triangle{3};
		
		U{i, j} = [triangle0(1), triangle1(1), triangle2(1)];
		V{i, j} = [triangle0(2), triangle1(2), triangle2(2)];
        W{i, j} = [triangle0(3), triangle1(3), triangle2(3)];
	end
end

X = cell2mat(X);
Y = cell2mat(Y);
Z = cell2mat(Z);
U = cell2mat(U);
V = cell2mat(V);
W = cell2mat(W);
end