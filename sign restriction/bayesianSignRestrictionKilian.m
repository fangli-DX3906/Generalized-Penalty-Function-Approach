function [IRF_q, theta_q, IRF_a, theta_a, IRF_c, theta_c] = bayesianSignRestrictionKilian(y, s, maxOrder, maxAuxOrder, nburn, nsim, M, H, signMAT, h1, h2)
       
    [T,N] = size(y);
    p = orderSelect(y, maxOrder);
    [A,~, ~,const,X,SIGMA] = olsvarc(y,p);	 % Frequentist approach: VAR with intercept
    Ns = size(s, 2);
    nShocks = size(signMAT,2);
    intval = nsim/10;
%     varD = 1: N;
%     varS = N+1:2*N;
%     auxD = 1:Ns;
%     auxS = Ns+1:2*Ns;  
%     [As, ps] = olsvaraux(s, y, 0, maxAuxOrder);

    % Check if the matrix algebra is correct
    Bhat        = [const(1:N,1), A(1:N,:)]'; 
    sigmautilde = SIGMA(1:N,1:N);
    Y           = y(1+p:end,:);
    X           = X';
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
    IRF_sub     = zeros(N*2*H, M*nsim);
%     IRF_q = zeros((N+Ns)*2,H, c);
%     IRF_a = zeros((N+Ns)*2,H, c);
    tot_draws  = M*nsim;

    for isim = 1:ntot

        % This is the covariance for the posterior density of alpha
        COV    = kron(sigmau,(X'*X)^(-1));

        % Posterior of beta|sigmau,Y,X ~ Normal(bethat,kron(sigmau,(X'X)^(1)))
        bet    = bethat + chol(COV)'*randn(N*(1+N*p),1);
        B      = (reshape(bet,1+N*p,N));

        % Draw from the posterior of sigmau|Y,X ~ IW(S,nu)
        sigmau = inv(wish(inv(S),nu));

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
                [Q,R] = qr(randn(N,N)); % draw orthogonal matrix Q
                Q = sign(diag(R)).*Q; % flip the sign of the columns of Q if elements along the main diagonal of R are negative
                candidate  = irfvar(BB,sigmau,p,H,Q); % candidate impulse response functions
                candsim = sum(candidate(:,h1:h2),2); % sign of the candidate impulse response functions. Also take care of the horizon
                candsimMAT = (reshape(candsim,N,N))';
                signsimMAT = sign(candsimMAT);
                if sum(sum((signsimMAT - signMAT).^2)) == 0
                    counter = counter + 1;
                    counter2 = counter2 + 1;
                    IRF_sub(:,counter2) = reshape(irfvar(BB,sigmau,p,H-1,Q), [], 1); 
                end
            end
        end 
    end
    
    loss_q = [];
    loss_a = [];
    loss_cosine = [];
    rs = N*2;
    for k=1:tot_draws
        theta_bar_mat = repmat(IRF_sub(:, k), 1, tot_draws);
        loss_q = [loss_q sum(sum((IRF_sub - theta_bar_mat).^2))/tot_draws];
        loss_a = [loss_a sum(sum(abs(IRF_sub - theta_bar_mat)))/tot_draws];
        
        % angular similarity
        theta_bar = reshape(IRF_sub(:, k), N*2, []);
        lc = 0;
        for h=1:tot_draws
            theta_i = reshape(IRF_sub(:, h), N*2, []);
            for r = 1: rs
                lc = lc + acos((theta_i(r,:)*theta_bar(r,:)')/(norm(theta_i(r,:))*norm(theta_bar(r,:))))/(rs^2*pi);
            end   
        end
        loss_cosine = [loss_cosine lc/tot_draws];
    end
    
    [~, idx_q] = min(loss_q);
    [~, idx_a] = min(loss_a);
    [~, idx_c] = min(real(loss_cosine));
    
    theta_q = reshape(IRF_sub(:, idx_q), N*2, H);
%     IRFaux_opt_q = recoverIRF(size(signMAT,2), N, Ns, theta_q, As, ps, H-1, 0);
%     IRFaux_opt_q = (reshape(IRFaux_opt_q, H, Ns*nShocks))';
%     IRF_q_opt=[theta_q(varD, :); IRFaux_opt_q(auxD,:); theta_q(varS,:); IRFaux_opt_q(auxS,:)];
    
    theta_a = reshape(IRF_sub(:, idx_a), N*2, H);
%     IRFaux_opt_a = recoverIRF(size(signMAT,2), N, Ns, theta_a, As, ps, H-1, 0);
%     IRFaux_opt_a = (reshape(IRFaux_opt_a, H, Ns*nShocks))';
%     IRF_a_opt=[theta_a(varD, :); IRFaux_opt_a(auxD,:); theta_a(varS,:); IRFaux_opt_a(auxS,:)];

    theta_c = reshape(IRF_sub(:, idx_c), N*2, H);
%     IRFaux_opt_c = recoverIRF(size(signMAT,2), N, Ns, theta_c, As, ps, H-1, 0);
%     IRFaux_opt_c = (reshape(IRFaux_opt_c, H, Ns*nShocks))';
%     IRF_c_opt=[theta_c(varD, :); IRFaux_opt_c(auxD,:); theta_c(varS,:); IRFaux_opt_c(auxS,:)];
    
    c = ceil(tot_draws*0.95);
    [~, cs_q] = sort(loss_q, 'ascend');
    [~, cs_a] = sort(loss_a, 'ascend');
    [~, cs_c] = sort(loss_cosine, 'ascend');
    cs_q = cs_q(1:c);
    cs_a = cs_a(1:c);
    cs_c = cs_c(1:c);

    for jj=1:c
        temp_q = reshape(IRF_sub(:, cs_q(jj)), N*2, H);
        temp_a = reshape(IRF_sub(:, cs_a(jj)), N*2, H);
        temp_c = reshape(IRF_sub(:, cs_c(jj)), N*2, H);

        IRF_q(:,:,jj) = temp_q;
        IRF_a(:,:,jj) = temp_a;
        IRF_c(:,:,jj) = temp_c;
        
%         IRFauxr_q = recoverIRF(size(signMAT,2), N, Ns, temp_q, As, ps, H-1, 0);
%         IRFauxr_q = (reshape(IRFauxr_q, H, Ns*nShocks))';            
%         IRF_q(:,:,jj)=[temp_q(varD, :); IRFauxr_q(auxD,:); temp_q(varS,:); IRFauxr_q(auxS,:)];
        
%         IRFauxr_a = recoverIRF(size(signMAT,2), N, Ns, temp_a, As, ps, H-1, 0);
%         IRFauxr_a = (reshape(IRFauxr_a, H, Ns*nShocks))';            
%         IRF_a(:,:,jj)=[temp_a(varD, :); IRFauxr_a(auxD,:); temp_a(varS,:); IRFauxr_a(auxS,:)];

%         IRFauxr_c = recoverIRF(size(signMAT,2), N, Ns, temp_c, As, ps, H-1, 0);
%         IRFauxr_c = (reshape(IRFauxr_c, H, Ns*nShocks))';            
%         IRF_c(:,:,jj)=[temp_c(varD, :); IRFauxr_c(auxD,:); temp_c(varS,:); IRFauxr_c(auxS,:)];
    end
   
end
