function hello = vector_chirality(infile, outfile, chifile)

%Calculate scalar chirality of a lattice of spins
%assume kagome lattice a1 = (1,0); a2 = (1/2, sqrt(3)/2)
%this means there is a basis - 3 spins per site, numbered
%counterclockwise
%        / \     / \     / 
%       10--11--13--14--16
%      /     \ /     \ /
%     3       6       9
%    / \     / \     / \
%   1---2---4---5---7---8
tic
raw_spins = csvread(infile);
spins = split_spins(raw_spins);%this is now a cell array of 3 3-vectors at each site
N = length(spins);%number of sites
L = sqrt(N);%length of each dimension
a1 = [1, 0];
a2 = [.5, sqrt(3)/2];
Sij = cell(L,L);
for j = 1:L
  for i = 1:L
    Sij{i,j} = spins{(i - 1)*L + j};
  end
end
chi = cell(2*(L - 1),(L - 1));
chi_sum = zeros(1, 3);
%This will be organized as follows: there are two indices for each lattice site
%They will each represent one triangle.  For example, indices 1 and 2 of the first
%column of chi will contain the chirality at A and at B, respectively:
%       10--11--13--14--16
%      /     \B/     \ /
%     3       6       9
%    /A\     / \     / \
%   1---2---4---5---7---8
for x = 1:(L - 1)
    for y = 1:(L - 1)
		%need to retrieve spins from 4 different triangles:
		triangle1 = Sij{x,y};%123
		triangle2 = Sij{x + 1, y};%456 - we need 6
		triangle3 = Sij{x, y + 1};%10,11,12 - we need 11
		triangle4 = Sij{x + 1, y + 1};%13,14,15 - we need 13
        %first triangle
        chi{2*x - 1, y} = cross(triangle1{1}, triangle1{2}) + cross(triangle1{2}, triangle1{3}) + cross(triangle1{3}, triangle1{1});
        %second triangle
        chi{2*x, y} = cross(triangle2{3}, triangle4{1}) + cross(triangle4{1}, triangle3{2}) + cross(triangle3{2}, triangle2{3});
        chi_sum = chi_sum + chi{2*x - 1, y} + chi{2*x, y};
    end
end

dims = size(chi);
avg_chi = chi_sum ./(dims(1)*dims(2));

chi_out = fopen(chifile, 'a');
fprintf(chi_out, '%f ,', avg_chi(1));
fprintf(chi_out, '%f ,', avg_chi(2));
fprintf(chi_out, '%f \n', avg_chi(3));
fclose(chi_out);

%create a quiver plot for the spins at each site
figure()
hold on
colormap(redblue)
[x,y,u,v] = quiver_kagome2D(a1,a2,Sij,L,L);
%quiver(x,y,u,v, .5);
%fill triangles
for i = 1:(L - 1)
	for j = 1:(L - 1)
		xtri1 = [x(i,3*j - 2),x(i,3*j - 1),x(i,3*j)];
		ytri1 = [y(i,3*j - 2),y(i,3*j - 1),y(i,3*j)];

		xtri2 = [x(i,3*j + 3),x(i + 1, 3*j - 1),x(i + 1, 3*j + 1)];%sites 6,11,13 as above
		ytri2 = [y(i,3*j + 3),y(i + 1, 3*j - 1),y(i + 1, 3*j + 1)];

		xhex = [x(i,3*j -1), x(i, 3*j), x(i + 1, 3*j - 2), x(i + 1, 3*j - 1), x(i, 3*j + 3), x(i, 3*j + 1)];%hexagon
		yhex = [y(i,3*j -1), y(i, 3*j), y(i + 1, 3*j - 2), y(i + 1, 3*j - 1), y(i, 3*j + 3), y(i, 3*j + 1)];%hexagon

		chitri1 = chi{2*i - 1, j};
		chitri2 = chi{2*i, j};
		if(i ~= 1 && i ~= (L - 1) && j ~= 1 && j ~= (L-1))
        chihex = (chitri1 + chitri2 + chi{2*i - 2, j} + chi{2*i + 1, j} + chi{2*i - 1, j + 1} + chi{2*i, j - 1}) ./ 6;
		else
			chihex = .5 * (chitri1 + chitri2);
		end

		fill(xtri1,ytri1,chitri1(3));
		fill(xtri2,ytri2,chitri2(3));
		fill(xhex, yhex, chihex(3));
	end
end
saveas(gcf,strcat(outfile,'.fig'))
saveas(gcf,strcat(outfile,'.png'))
hold off

toc
end