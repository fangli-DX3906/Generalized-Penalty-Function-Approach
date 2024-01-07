/* Feb 11th, 2022 */
/* NK with capital */

var Y C N w I m lam pi A T chi K r betav;
varexo e_a e_t e_chi e_beta;
parameters alpha bet delta sigmac sigman delta eta pp gammap rhoa rhot rhoi rhoc rhob psi;

alpha = 0.33;
bet = 0.99;
delta = 0.025;
sigmac = 2;
sigman = 2;
eta = 2;
pp = 1.5;
rhoa = 0.95;
rhot = 0.9;
rhoi = 0.75;
rhoc = 0.9;
rhob = 0.75;
gammap = 1;

% calibrate psi
[SS, psi] = NKSS(alpha, bet, sigmac, sigman, delta, eta, 0.33, 0.33*ones(1,9));

model;
    [name = 'intratemporal']
    w = C^sigmac * psi * (1-N)^sigman;
    
    [name = 'Euler']
    1 = betv * (C(+1)/C)^(-sigmac) * (1+r)/pi(+1);
%     1 = bet * (C(+1)/C)^(-sigmac) * (1+r(-1))/pi(+1);

    [name = 'production']
    Y = A * K(-1)^alpha * L^(1-alpha);
    
    [name = 'TFP shock']
    log(A) = rhoa * log(A(-1)) + e_a;
    
    [name = 'FOC P']
    1 - gammap * (pi - SS(8)) * pi - eta*lam_t + m*gammap*(pi(+1)-SS(8))*pi(+1)*Y(+1) / Y = 0; 
    
    [name = 'FOC K']
    1 = m*(1-lam(+1))*alpha*A(+1)*K^(alpha-1)*N(+1)^(1-alpha)+1-delta;
    
    [name = 'FOC N']
    (1-lam)*(1-alpha) = A*K(-1)^(alpha)*N^(-alpha);
    
    [name = 'Monetary policy']
    (1+r(-1)) = (1+r(-2))^rhoi * ((1+SS(13))*(pi/SS(8))^pp)^(1-rhoi)*chi;
    
    [name = 'Monetary policy shock']
    log(chi) = rhoc * log(chi(-1)) + e_chi;
    
    [name = 'LMO']
    K = (1-delta)*K(-1) + I;
    
    [name = 'market clearing']
    Y = C + T + I + 0.5*gammap*(pi-SS(8))^2 * Y;
    
    [name = 'spending shock']
    T = rhot * T(-1) + e_t;
    
    [name = 'SDF']
    m = betv * (C(+1)/C)^sigmac;
    
    [name = 'beta shock']
    log(betv/bet) = rhob*log(betv(-1)/bet) + e_beta;
end;

steady_state_model;
    Y = SS(1);
    C = SS(2);
    N = SS(3);
    w = SS(4);
    I = SS(5);
    m = SS(6);
    lam = SS(7);
    pi = SS(8);
    A = SS(9);
    T = SS(10);
    chi = SS(11);
    K = SS(12);
    r = SS(13);
    betv = SS(14);
end;
steady;
check;
resid;

shocks;
    var e_a = 0.01;
    var e_t = 0.01;
    var e_chi = 0.01;
    var e_beta = 0.01
end;
    
stoch_simul(periods = 1000, irf = 40);