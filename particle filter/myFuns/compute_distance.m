function distance = compute_distance(x0,y0,phi,contour)
    % define frame for alpha and beta
    alpha = zeros(1,10);
    beta  = zeros(1,10);
    contour = [contour;
               contour(1,:)];
    % read corner coordinate
    for i = 1:10
        x1 = contour( i ,1);
        y1 = contour( i ,2);
        x2 = contour(i+1,1);
        y2 = contour(i+1,2);
        det = sin(phi)*(x2-x1)+cos(phi)*(y1-y2);
        alpha(i) = (sin(phi)*(x2-x0)+cos(phi)*(y0-y2))/det;
        beta(i)  = ( (y2-y1)*(x0-x2)+(x1-x2)*(y0-y2) )/det;
    end
    % determine the wall intersected
    W = alpha >= 0 & alpha <= 1 & beta >= 0;
    B = beta.*W;
    distance = min(B(B~=0));
    if isempty(distance)
        distance = inf;
    end
end