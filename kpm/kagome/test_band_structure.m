%Plot a test band structure
k_range = linspace(4*pi/3, 0, 1000);

energies = zeros(length(k_range), 24);
i = 1;
for k = k_range
    k_in = [k, 0];
    ham = band_structure(k_in);
    energies(i, :) = eig(ham);
    i = i + 1;
end
figure()
hold on
for i = 1:24
    plot(energies(:, i))
end
