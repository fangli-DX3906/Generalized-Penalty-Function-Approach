function IRF = bayesianSignRestriction(y, s, maxOrder, maxAuxOrder, nburn, nsim, M, H, signMAT, h1, h2)
       
    [T,N] = size(y);
    [p, p_1, p_2] = selectOptimLagOrder(y, maxOrder);   % for now
%     p = orderSelect(y, maxOrder);
    [A,~, ~,const,X,SIGMA] = olsvarc(y,p);	 % Frequentist approach: VAR with intercept
    Ns = size(s, 2);
    nShocks = size(signMAT,2);
    intval = nsim/10;

    % Check if the matrix algebra is correct
    Bhat        = [const(1:N,1), A(1:N,:)]'; % re-add constant and remove identity matrix from the companion form
    sigmautilde = SIGMA(1:N,1:N);
    Y           = y(1+p:end,:);
    X           = X';

    % Check if the vectorization algebra is correct
%     ybar     = reshape(y(1+p:end,:),[],1);
%     ubar     = reshape(Uhat(1:N,:)',[],1);
    bethat   = reshape(Bhat,[],1);

    % Posterior distribution parameters
    S        = T*sigmautilde;  % Inverse Wishart parameter 1
    if p > 1
        nu       = (T-N)*(p-1);    % Inverse Wishart parameter 2
    else
        nu       = T-N;    % Inverse Wishart parameter 2
    end

    % Simulate impulse response functions
    ntot           = nsim + nburn;   
    sigmau         = sigmautilde; % initial value for sigmau. Then it will use the drawns
    counter2       = 0;

    % Pre-allocate matrices
    IRF            = zeros((N+Ns)*2,H,M*nsim);
    % Not calculating the variance decomposition and historical decomposition
%     VD             = zeros(N^2,H+1,M*nsim);
%     vd             = zeros(N,N,H+1);
%     yhat           = zeros(N,T-p,N,M*nsim); 

    for isim = 1:ntot

%         disp(['Simulation ',num2str(isim)])

        % This is the covariance for the posterior density of alpha
        COV    = kron(sigmau,(X'*X)^(-1));

        % Posterior of beta|sigmau,Y,X ~ Normal(bethat,kron(sigmau,(X'X)^(1)))
        bet    = bethat + chol(COV)'*randn(N*(1+N*p),1);
        B      = (reshape(bet,1+N*p,N));
%         usim   = Y - X*B;

        % Draw from the posterior of sigmau|Y,X ~ IW(S,nu)
        sigmau = inv(wish(inv(S),nu));
        dsigma = inverseWishartDist(S, nu); 


        if isim > nburn
            if mod((isim-nburn),intval)==0
                completePct = 100*(isim-nburn)/nsim;
                disp(['   ', num2str(completePct), '% completed!'])
            end
            
            % Calculate IRFs
            id = eye(N*p);
            BB = [B(2:end,:)'; id(1:end-N,:)];
            counter = 0;
            while counter < M 
                [Q,R]               = qr(randn(N,N)); % draw orthogonal matrix Q
                Q                   = sign(diag(R)).*Q; % flip the sign of the columns of Q if elements along the main diagonal of R are negative
                candidate           = irfvar(BB,sigmau,p,H,Q); % candidate impulse response functions
                candsim             = sum(candidate(:,h1:h2),2); % sign of the candidate impulse response functions. Also take care of the horizon
                candsimMAT          = (reshape(candsim,N,N))';
                signsimMAT          = sign(candsimMAT);
                if sum(sum((signsimMAT - signMAT).^2)) == 0
                    counter           = counter + 1;
                    counter2          = counter2 + 1;

                    % Calculate impulse responses
                    IRFr = irfvar(BB,sigmau,p,H-1,Q); 
                    [As, ps] = olsvaraux(s, y, 0, maxAuxOrder);
                    IRFauxr = recoverIRF(size(signMAT,2), N, Ns, IRFr, As, ps, H-1, 0);
                    IRFauxr = (reshape(IRFauxr, H, Ns*nShocks))';
                    varD = 1: N;
                    varS = N+1:2*N;
                    auxD = 1:Ns;
                    auxS = Ns+1:2*Ns;               
                    IRF(:,:,counter2)=[IRFr(varD, :); IRFauxr(auxD,:); IRFr(varS,:); IRFauxr(auxS,:)];

                    % Calculate variance decomposition
%                     IRsq              = reshape(IRF(:,1:H+1,counter2).^2,N,N,H+1);
%                     IRsqsumH          = cumsum(IRsq,3);
%                     totFEV            = squeeze(sum(IRsqsumH,2));
%                     for ih = 1:size(totFEV,2)
%                         vd(:,:,ih) = IRsqsumH(:,:,ih)./totFEV(:,ih);
%                     end
%                     VD(:,:,counter2)  = reshape(vd,N^2,H+1);

%                     % Calculate historical decomposition
%                     shocksim          = (chol(sigmau)'*Q) \ (usim'); %compute structural shocks shocksim from reduced form shocks usim
%                     for it = 1:T-p
%                         for is = 1:N
%                             for iv = 1:N
%                                 yhat(iv,it,is,counter2) = dot(IRF(iv+(is-1)*N,1:it,counter2),shocksim(is,it:-1:1));
%                             end
%                         end
%                     end

                end
            end
        end 
    end

end
