function [x_r,y_r,phi] = random_sampling(N)
    P = rand(3,2*N);
    P(1,:) = P(1,:) * 2;
    P(2,:) = P(2,:) * 2;
    P(3,:) = P(3,:) * 2 * pi;
    for i = 1:2*N
        x = P(1,i);
        y = P(2,i);
        if (x<=0.5&&y<=1.0) || (x<=0.25&&y>=1.75) || (x>=1.5&&y>=1.0) || (x>=1.5&&y>=-x+2.5)
            P(:,i) = [0;0;0];
        end
    end
    P = nonzeros(P);
    P = reshape(P',3,[]);
    x_r = P(1,1:N);
    y_r = P(2,1:N);
    phi = P(3,1:N);
end