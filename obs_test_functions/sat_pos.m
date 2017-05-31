function [rk,dt_rel] = sat_pos( gps_ephemeris, Secs)

mu = 3.986005e14; % Gravitational parameter [m^/s^]
omega_Earth = 7.2921151467e-5;% Earth's roration [rad/s]
c = 299792458; % speed of light [m/s]

% Extract all ephemeris components
ephem = num2cell(gps_ephemeris);
size(ephem);
%ephem = gps_ephemeris;
[prn,M0,delta_n,ecc,sqrt_a,Loa,incl,perigee,ra_rate,i_rate,Cuc,Cus,Crc,Crs,Cic,Cis,...
    Toc,IODE,GPS_week,Toc,Af0,Af1,Af2,nil,health] = deal(ephem{:});
A = sqrt_a^2;
% Correct Mean Motion
n0  = sqrt(mu/(A)^3); % Calculated mean motion [rad/s]
n   = n0 + delta_n;                 % Corrected Mean Motion

% Correct Time
tk  = Secs-Toc;  

% Mean Anomaly
Mk = M0 + n*tk; % Mean anomaly

%------Eccentric Anomaly
options=optimset('Display','off','TolFun',1e-15,'TolX',1e-15);
Ek = fsolve(@(Ek) (Ek)-ecc*sin(Ek)-Mk,4,options);

% True Anomaly
vk = atan2(      (sqrt(1-ecc^2)*sin(Ek)/(1-ecc*cos(Ek))), ...
                      ((cos(Ek)-ecc)/(1-ecc*cos(Ek)))   );
                  
% Argument of Latitude
Phik = vk + perigee;

% Second Harmonic Perturbations
del_uk = Cus*sin(2*Phik) + Cuc*cos(2*Phik);
del_rk = Crs*sin(2*Phik) + Crc*cos(2*Phik);
del_ik = Cis*sin(2*Phik) + Cic*cos(2*Phik);

% Corrected argumet of latitude, radius, inclination
uk = Phik + del_uk;
rk = A*(1-ecc*cos(Ek)) + del_rk;
ik = incl + del_ik + i_rate*tk;

% Position in Orbit Plane
xkp = rk*cos(uk);
ykp = rk*sin(uk);

% Corrected Longitude of ascending node
Omegak = Loa + (ra_rate - omega_Earth)*tk - omega_Earth*Toc;

% Earth Fixed Coordinates
xk = xkp * cos(Omegak) - ykp * cos(ik) * sin(Omegak);
yk = xkp * sin(Omegak) + ykp * cos(ik) * cos(Omegak);
zk = ykp * sin(ik);

% Relativity time shift
dt_rel = 2*sqrt(mu)/c^2 * ecc * sqrt_a * sin(Ek);
rk = [xk,yk,zk];

