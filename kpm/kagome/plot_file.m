function data = plot_file(infile)
%simple function to plot data from a file
data = csvread(infile);
plot(data(1:(end - 1)))

end