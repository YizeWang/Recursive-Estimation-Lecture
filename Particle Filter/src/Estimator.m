function postParticles = Estimator(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==1, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the 
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                           
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index, scalar
%                       corresponds to continous time t = k*Ts
%                       If tm==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]

N = 1000; % number of particles
M = 10; % massive sampling coefficient
K = 0.1; % roughening parameter
WarningFlag = 1; % whether display meansurement update failure

%% Initialization
if (km == 0)
    start_region = rand(1,N) < 0.5; % 0 stands for A, 1 stands for B
    r = sqrt( rand(1,N) ) * estConst.d; % 1xN matrix, random polar distance
    theta = rand(1,N) * 2 * pi; % 1xN matrix, random polar angle
    x_relative = r .* cos(theta); % 1xN matrix, relative x coordinate to the center
    y_relative = r .* sin(theta); % 1xN matrix, relative y coordinate to the center
    postParticles.x_r = x_relative + ( start_region * (estConst.pB(1)-estConst.pA(1)) + estConst.pA(1) ); % 1xN matrix, particle x coordinate
    postParticles.y_r = y_relative + ( start_region * (estConst.pB(2)-estConst.pA(2)) + estConst.pA(2) ); % 1xN matrix, particle y coordinate
    postParticles.phi = rand(1,N) * 2 * estConst.phi_0 - estConst.phi_0; % 1xN_particles matrix
    return;
end

%% Prior Update
% prior update particles through process model
vf   = ( rand(1,N) - 0.5 ) * estConst.sigma_f;
vphi = ( rand(1,N) - 0.5 ) * estConst.sigma_phi;
priorParticles.x_r = prevPostParticles.x_r + (act(1)+vf).*cos(prevPostParticles.phi);
priorParticles.y_r = prevPostParticles.y_r + (act(1)+vf).*sin(prevPostParticles.phi);
priorParticles.phi = prevPostParticles.phi + act(2) + vphi;

%% Posterior Update
% store prior update distance
distance = zeros(1,N);
beta = zeros(1,N);
for k = 1:N
    distance(k) = compute_distance(priorParticles.x_r(k),priorParticles.y_r(k),priorParticles.phi(k),estConst.contour);
    beta(k) = read_noise_probability(sens-distance(k),estConst.epsilon);
end
% particle weights
Normalization = sum(beta);
while(~Normalization)
    if WarningFlag
        string = "Zero Probability after Measurement Update at " + int2str(km) + ". Resampling ...";
        disp(string);
    end
    [priorParticles.x_r,priorParticles.y_r,priorParticles.phi] = random_sampling(M*N); % massive sampling
    distance = zeros(1,N);
    beta = zeros(1,N);
    for k = 1:M*N
        distance(k) = compute_distance(priorParticles.x_r(k),priorParticles.y_r(k),priorParticles.phi(k),estConst.contour);
        beta(k) = read_noise_probability(sens-distance(k),estConst.epsilon);
    end
    [~,Ind] = maxk(beta,N); % find N most probable particles
    % truncate particles
    priorParticles.x_r = priorParticles.x_r(Ind);
    priorParticles.y_r = priorParticles.y_r(Ind);
    priorParticles.phi = priorParticles.phi(Ind);
    beta = beta(Ind);
    Normalization = sum(beta);
end
beta = beta / Normalization;
% resampling
bin = resampling(beta);
for k = 1:N
    postParticles.x_r(k) = priorParticles.x_r(bin(k));
    postParticles.y_r(k) = priorParticles.y_r(bin(k));
    postParticles.phi(k) = priorParticles.phi(bin(k));
end
% roughening
[postParticles.x_r,postParticles.y_r,postParticles.phi] = roughening(K,N,postParticles.x_r,postParticles.y_r,postParticles.phi);
end