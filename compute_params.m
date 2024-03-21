function [mu, SIGMA] = compute_params(t, alph, noise_var_linear, PHI)
  
    A = diag(alph); %diagonal matrix of alphas'.

    SIGMA = ((1/noise_var_linear).*PHI'*PHI + A)^(-1);
    mu = (1/noise_var_linear)*SIGMA*PHI'*t;
 


    for j = 1:500

        muold = mu; %to compare previous value with current value.
    
        %calculating alpha_i's
        for i = 1:length(alph)
            gammai = 1 - alph(i)*SIGMA(i,i);
            alph(i) = gammai/(mu(i)^2);
            
        end
    
        %updating values of SIGMA and mu.
        A = diag(alph);
        SIGMA = ((1/noise_var_linear)*(PHI'*PHI) + A)^(-1);
        mu = (1/noise_var_linear)*SIGMA*PHI'*t;
    
        %breaking condition
        if((norm(mu - muold)/norm(muold))^2 <= 0.001)
            break;
        end
      
    end

end