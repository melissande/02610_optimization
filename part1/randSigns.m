function K = randSigns(x)

%Remove all zero entries
x(x==0) = [];

%Count n+, n and m
m = length(x);
n_p = sum(x>0);
n_m = sum(x<0);

%Calculate the mean and standard deviation of the runs

mu = (2*n_p*n_m / m) +1;
s = sqrt ( (mu-1)*(mu-2) / (m-1) );

%Create a boolean
b = x > 0;

%Count the number of runs in the boolean (i.e. in x)
N = 0;
    for i = 1:length(b)-1
        if b(i) ~= b(i+1)
            N = N+1;
        end  
    end

Z = norm(N-mu) / s;

if Z <= 1.96
    fprintf('The residual sign sequence appears random with Z = %f, which is less than 1.96',Z);
    K = Z;
else
    fprintf('The residual sign sequence does NOT appear random with Z = %f, which is larger than 1.96',Z);
    K = Z;
end

end




