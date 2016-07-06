function s = split_spins(spins)
%Spins is the raw spin data for {s1, s2, s3...sN}; want to 
%split it from a 3N vector to an N vector of 3-vectors
%assume length(spins) is divisible by 3
N = length(spins)/3;
s = cell(1, N);
for i = 1:N
	s{i} = [spins(3*i - 2), spins(3*i - 1), spins(3*i)];
end
end
