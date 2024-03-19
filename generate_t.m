function [t, w, PHI, epsilon] = generate_t(M, N, D0, noise_var_linear)

    t = zeros(N,1);

    PHI = normrnd(0,1,N,M); %generating PHI/design matrix.

    %generating weights
    w = zeros(M,1);
    rand_indices = randperm(M,D0); %generating random indices.
    for i = rand_indices
        w(i) = normrnd(0,sqrt(noise_var_linear));
    end
    % disp(w);
    
    epsilon = normrnd(0,sqrt(noise_var_linear),N,1); %noise samples.
    
    %generating t
    t = PHI*w + epsilon;


end