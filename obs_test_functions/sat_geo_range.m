
function [R, r1, r2 ,r3 , dt_rel] = sat_geo_range(rStation, GPS_Weeks, GPS_SOW, PRN, epoch_nav_data, Secs)
mu = 3.986005e14; % Gravitational parameter [m^/s^]
omega_Earth = 7.2921151467e-5;% Earth's roration [rad/s]
c = 299792458; % speed of light [m/s]


% Get Sat Position from Ephemeris data
[rSat,dt_rel] = sat_pos( epoch_nav_data, Secs);


% Set up convergence limits
R = 0;
conv_limit = 1e-12;
max_iters = 100;
iter = 1;

% Iterate and converge on Geometric Range
while(1)
    % Calculate Geometric Range
    Rtmp = norm( rSat - rStation );
    
    % Check for Convergence
    if(abs(Rtmp - R) < conv_limit)
        break
    end
    
    % Assign new Range Value now that criterion are passed
    R = Rtmp;
    
    % Check for iteration limit
    if(iter > max_iters)
        error('Range Calculation not converging!')
    end
    
    % Increase iteration count
    iter = iter + 1;
    
    % Calculate 'Tt', time of transmission
    dt = R/c;
%   Tr = epochData(SOW_column);
    Tr = GPS_SOW;
    Tt = Tr - dt;
    
    % Recalculate Satellite position
    Secs = Tt;     
    % to use new time value
    [rSat,dt_rel] = sat_pos(epoch_nav_data,Secs);
    
    % Rotate Sat position at time Tr (account for earth's rotation)
    phi = omega_Earth*dt;
    rot3 = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
    rSat = transpose(rot3*rSat');
    
    % Calculate the new range
    
    R = norm(rSat - rStation);
    r1 = rSat(1);
    r2 = rSat(2);
    r3 = rSat(3);
end
