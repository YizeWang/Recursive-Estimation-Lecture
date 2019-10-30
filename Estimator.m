function [posEst, linVelEst, oriEst, driftEst, ...
          posVar, linVelVar, oriVar, driftVar, estState] = ...
          Estimator(estState, actuate, sense, tm, estConst)
% [posEst, linVelEst, oriEst, driftEst,
%  posVar, linVelVar, oriVar, driftVar, estState] = 
%  Estimator(estState, actuate, sense, tm, estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.

%% Initialization
if (tm == 0)
    % initialization
    
    % initial state mean
    posEst = [0 0]; % 1x2 matrix
    linVelEst = [0 0]; % 1x2 matrix
    oriEst = 0; % 1x1 matrix
    driftEst = 0; % 1x1 matrix
    
    % initial state variance
    posVar = [estConst.StartRadiusBound^2/4 estConst.StartRadiusBound^2/4]; % 1x2 matrix
    linVelVar = [0 0]; % 1x2 matrix
    oriVar = estConst.RotationStartBound^2/3; % 1x1 matrix
    driftVar = 0; % 1x1 matrix
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar oriVar linVelVar driftVar]);
    % estimator state init (initial posterior state)
    estState.xm = [posEst oriEst linVelEst driftEst]; % 6*1 matrix
    % time of last update
    estState.tm = tm;
    return;
end
%% Estimator iteration.

% get time since last estimator update
dt = tm - estState.tm;

% update measurement time
estState.tm = tm;

% define constant matrix Qc and R
Qc = diag([estConst.DragNoise estConst.RudderNoise estConst.GyroDriftNoise]);
R  = diag([estConst.DistNoiseA estConst.DistNoiseB estConst.DistNoiseC estConst.GyroNoise estConst.CompassNoise]);

% define A and L matrix
A = @(x)  [0 0 0 1 0 0;
           0 0 0 0 1 0;
           0 0 0 0 0 0;
           0 0 -sin(x(3)) * ( tanh(actuate(1)) - estConst.dragCoefficient*( x(4)^2 + x(5)^2 ) ) -2*estConst.dragCoefficient*x(4)*cos(x(3)) -2*estConst.dragCoefficient*x(5)*cos(x(3)) 0;
           0 0  cos(x(3)) * ( tanh(actuate(1)) - estConst.dragCoefficient*( x(4)^2 + x(5)^2 ) ) -2*estConst.dragCoefficient*x(4)*sin(x(3)) -2*estConst.dragCoefficient*x(5)*sin(x(3)) 0;
           0 0 0 0 0 0];
L = @(x)  [0 0 0;
           0 0 0;
           0 estConst.rudderCoefficient * actuate(2) 0;
           -cos(x(3)) * estConst.dragCoefficient * ( x(4)^2 + x(5)^2 ) 0 0;
           -sin(x(3)) * estConst.dragCoefficient * ( x(4)^2 + x(5)^2 ) 0 0;
           0 0 1];
       
% reshape P matrix
P = @(x) reshape(x(7:end),[6,6]); % 6*6 matrix
Pdot = @(x) reshape(A(x)*P(x) + P(x)*A(x)' + L(x)*Qc*L(x)',[36,1]);
xdot = @(x) [x(4);
             x(5);
             estConst.rudderCoefficient * actuate(2);
             cos(x(3))*( tanh(actuate(1)) - estConst.dragCoefficient*( x(4)^2 + x(5)^2) );
             sin(x(3))*( tanh(actuate(1)) - estConst.dragCoefficient*( x(4)^2 + x(5)^2) );
             0];

% prior update

% solve differential equation
dxPdt = @(t,xP) [xdot(xP);
                 Pdot(xP)];
[~, xP] = ode45(dxPdt, [tm-dt tm], [estState.xm'; reshape(estState.Pm,[36,1])]);
xp = xP(end, 1:6)';
Pp = reshape( xP(end,7:42), [6,6] );

% measurement update

% define M matrix
M = eye(5);

% define h(\hat{x}_p,0)
da = @(x) sqrt( ( x(1) - estConst.pos_radioA(1) )^2 + ( x(2) - estConst.pos_radioA(2) )^2 );
db = @(x) sqrt( ( x(1) - estConst.pos_radioB(1) )^2 + ( x(2) - estConst.pos_radioB(2) )^2 );
dc = @(x) sqrt( ( x(1) - estConst.pos_radioC(1) )^2 + ( x(2) - estConst.pos_radioC(2) )^2 );
h  = @(x) [da(x);
           db(x);
           dc(x);
           x(3) + x(6);
           x(3)];
h = h(xp);

% define H matrix
H11 = @(x) ( x(1) - estConst.pos_radioA(1) ) / da(x);
H12 = @(x) ( x(2) - estConst.pos_radioA(2) ) / da(x);
H21 = @(x) ( x(1) - estConst.pos_radioB(1) ) / db(x);
H22 = @(x) ( x(2) - estConst.pos_radioB(2) ) / db(x);
H31 = @(x) ( x(1) - estConst.pos_radioC(1) ) / dc(x);
H32 = @(x) ( x(2) - estConst.pos_radioC(2) ) / dc(x);
H   = @(x) [H11(x) H12(x) 0 0 0 0;
            H21(x) H22(x) 0 0 0 0;
            H31(x) H32(x) 0 0 0 0;
                 0      0 1 0 0 1;
                 0      0 1 0 0 0];
H = H(xp);

% measurement is available at station c
if ~isinf(sense(3))
    disp(['measurement available at ' num2str(tm)]);
% measurement is not available at station c
else
    H(3,:) = [];
    M(3,:) = [];
    M(:,3) = [];
    R(3,:) = [];
    R(:,3) = [];
    sense(3) = [];
    h(3) = [];
end

% measurement update
K  = Pp * H' * inv( H * Pp * H' + M * R * M' );
xm = xp + K * (sense' - h);
Pm = ( eye(6) - K * H ) * Pp;

% Get resulting estimates and variances
% Output quantities
estState.xm = xm';
estState.Pm = Pm;
posEst = xm(1:2)';
linVelEst = xm(4:5)';
oriEst = xm(3);
driftEst = xm(6);

posVar = [Pm(1,1) Pm(2,2)];
linVelVar = [Pm(4,4) Pm(5,5)];
oriVar = Pm(3,3);
driftVar = Pm(6,6);

end