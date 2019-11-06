function probability = read_noise_probability(w,epsilon)
    w = abs(w);
    if w <= 2*epsilon
        probability =   -w/5/epsilon^2 + 2/5/epsilon;
    elseif w <= 2.5*epsilon
        probability =  2*w/5/epsilon^2 - 4/5/epsilon;
    elseif w<= 3*epsilon
        probability = -2*w/5/epsilon^2 + 6/5/epsilon;
    else
        probability = 0;
    end
end