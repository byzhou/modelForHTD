% This is some basic definitions, the size of the test matrix
dimension       = 5 ;
noiseStrength   = 0.1 ;
signalAve       = 1 ;
signalVar       = 0.5 ;

% This part is for the generation of the matrices
x_noiseFree     = normrnd (signalAve, signalVar, dimension, dimension);

% Real Monte Carlo Simulation - This is for the reference of the critiria
% distribution

testNum         = 1000 ;
% For the first rho
x_noise         = normrnd (0, noiseStrength, dimension, dimension);
x_withNoise     = x_noiseFree + x_noise ;
covariance      = cov (x_noiseFree , x_withNoise);
rho             = covariance(1 , 2) / sqrt (covariance(1 , 1) * covariance(2 , 2)) ;

% I use the end+1 to add more rho values into the rho queue
for i = 1 : (testNum - 1)
    x_noise         = normrnd (0, noiseStrength, dimension, dimension);
    x_withNoise     = x_noiseFree + x_noise ;
    covariance      = cov (x_noiseFree , x_withNoise);
    rho(end + 1)    = covariance(1 , 2) / sqrt (covariance(1 , 1) * covariance(2 , 2)) ; 
end

% show the histogram of the distribution
hist(rho , [0.9:0.01:1]);
