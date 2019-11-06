function [x_r,y_r,phi] = roughening(K,N,x_r,y_r,phi)
    E1 = max(x_r) - min(x_r);
    E2 = max(y_r) - min(y_r);
    E3 = max(phi) - min(phi);
    sigma1 = K * E1 / N^(1/3);
    sigma2 = K * E2 / N^(1/3);
    sigma3 = K * E3 / N^(1/3);
    x_r = x_r + sigma1*randn(1,N);
    y_r = y_r + sigma2*randn(1,N);
    phi = phi + sigma3*randn(1,N);
end

