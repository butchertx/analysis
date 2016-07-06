function s = split_spins(spins)
%Spins is the raw spin data for {s1, s2, s3...sN}; want to 
%split it from a 9N vector to an N vector of 3 3-vectors each
%assume length(spins) is divisible by 9
N = length(spins)/9;
s = cell(1, N);
for i = 1:N
	s{i} = {[spins(9*i - 8), spins(9*i - 7), spins(9*i - 6)],...
		[spins(9*i - 5), spins(9*i - 4), spins(9*i - 3)],...
		[spins(9*i - 2), spins(9*i - 1), spins(9*i)]};
end
end
