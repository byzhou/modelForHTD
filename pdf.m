% The normal distribution of alpha, E[RX] distribution
syms sigma_normal ...   % Standard deviation in normal distribution
     mu_normal ...      % Expectation in normal distribution
     pi ...             % Constant pi
     pdfnorm ...        % Normal distribution pdf
 ... % E[RX] distribution related variables
     x ...              % \rho value in FALSE POSITIVE
     sigma_x ...        % Standard deviation in \rho distribution
     mu_x ...           % Expectation in \rho distribution
     n ...              % Number of Monte-Carlo Simulation
     alpha ...          % random variable here, representing E[RX]
 ... % sigma_r distribution related variables
     beta ...           % random variable here, representing \sigma_{R}^2
 ... % constructed varaibles
     V ...              % actually used
     W ...              % not useful
 ... % jacobian determinant
     jacobi

mu_normal       = mu_x ;
sigma_normal    = 1 / n * ( sigma_x ^ 2 + mu_x ^ 2 ) ;

% alpha distribution
pdfnorm         = 1 / (sigma_normal * sqrt ( 2 * pi ) ) ...
                    * exp ( - (alpha - mu_normal) ^ 2 / ( 2 * sigma_normal ^ 2 ) * sqrt (2) ) ;

% beta distribution
pdfchisquare    = 1 / (2 ^ n * factorial (n - 1)) * beta ^ (n - 1) * exp ( - beta / 2 ) ;

% construct non-existing variable
% V               = ( sigma_x ^ 2 + alpha ) ^ 2 ...
%                     / ( sigma_x ^ 2 * ( sigma_x ^ 2 + beta * 2 + 2 * alpha ) ) ;
% W               = sigma_x ^ 2 * ( sigma_x ^ 2 + beta * 2 + 2 * alpha ) ;

% reverse representation
alpha           = sqrt ( V * W ) - sigma_x ^ 2 ;
beta            = W / sigma_x ^ 2 + sigma_x ^ 2 - 2 * sqrt ( V * W ) ;

% Jacobian determinant
jacobi          = diff ( alpha , V ) * diff ( beta , W ) - diff ( alpha , W ) * diff ( beta , V ) ;

% final distribution
mypdf           = pdfnorm * pdfchisquare * jacobi 
