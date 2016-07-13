function hello = scalar_chirality(infile, outfile, chifile)

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
chi = zeros(2*(L - 1),(L - 1));
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
        chi(2*x - 1, y) = dot(triangle1{1}, cross(triangle1{2}, triangle1{3}));
        %second triangle
        chi(2*x, y) = dot(triangle2{3}, cross(triangle4{1}, triangle3{2}));
  	end
end
%the points will be xval, yval, chirality, where xval,yval gives the leftmost corner of the triangle
%(for example, the first two points will be the locations of site 1 and site 11)
plot_mat = zeros(2 * (L - 1) * (L - 1), 3);
for i = 1:(L-1)
  for j = 1:(L-1)
    plot_mat(2*(i-1)*(L-1) + 2*j - 1, 1:2) = j*a1 + i*a2;%the first triangle
    plot_mat(2*(i-1)*(L-1) + 2*j - 1, 3) = chi(2*i - 1,j);
    plot_mat(2*(i-1)*(L-1) + 2*j, 1:2) = (j + 0.5)*a1 + (i + 1)*a2;%the second triangle
    plot_mat(2*(i-1)*(L-1) + 2*j, 3) = chi(2*i, j); 
  end
end
dims = size(chi);
avg_chi = sum(sum(chi))/(dims(1)*dims(2));
chi_max = max(max(abs(chi)));
%chi = chi / chi_max;%now chi takes values between -1 and 1 (for color mapping)

chi_out = fopen(chifile, 'a');
fprintf(chi_out, '%f , ', avg_chi);
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

		chitri1 = chi(2*i - 1, j);
		chitri2 = chi(2*i, j);
		if(i ~= 1 && i ~= (L - 1) && j ~= 1 && j ~= (L-1))
			chihex = (chitri1 + chitri2 + chi(2*i - 2, j) + chi(2*i + 1, j) + chi(2*i - 1, j + 1) + chi(2*i, j - 1))/6;
		else
			chihex = .5 * (chitri1 + chitri2);
		end

		fill(xtri1,ytri1,chitri1);
		fill(xtri2,ytri2,chitri2);
		fill(xhex, yhex, chihex);
	end
end
saveas(gcf,strcat(outfile,'.fig'))
saveas(gcf,strcat(outfile,'.png'))
hold off

toc
end
