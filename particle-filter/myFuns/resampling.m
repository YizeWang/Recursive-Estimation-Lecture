function bin = resampling(beta)
    beta = [0 beta];
    N = size(beta,2);
    gamma = cumsum(beta);
    [~,~,bin] = histcounts(rand(1,N),gamma);
end

