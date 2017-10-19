%Investigation for which N in the model that yields white noise errors under random signs
%test
stat.A = [];
stat.B = [];
stat.N = [];
stat.Z = [];

%Data
Y = [91.2200 28.0400 22.9100 26.6500 42.9600 101.0500 202.3600 328.0200 364.1200 299.2300 238.0000 227.4900 218.0300 223.6200 238.7500 271.2600 267.7200 251.3200 230.0400 206.6900 170.7700 131.6700 143.8500 157.5700];
T = 1:24;

%Residuals from NOfit are saved in stat.B, stat.N is the corresponding
%value of N while stat.Z is the random signs test for the particular
%residuals
for i=2:100
    N = 2*i-1;
    [A, B] = NOfit_Phil(T,Y,N);
     stat.B = [stat.B, B];
     stat.N = [stat.N, N];
     stat.Z = [stat.Z , randSigns(stat.B(:,i-1))];
end
%Plot Z-value as a function of N, with the line Z = 1.96
%Z values under this line confirms white noise data errors 
K = ones(1,length(stat.N))*1.96;
figure(2)
plot(stat.N,stat.Z)
hold on
plot(stat.N,K)