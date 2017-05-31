clear all; close all; clc

addpath obs
addpath obs_test_functions

observation_data = [];
ephemeris_data = [];
ionospheric_data = [];

% Setting a high enough quantity for LS iteration
dx=100; dy=dx; dz=dx;

% Gravitaional paramerter, rotation rate and speed of light
mu = 3.986005e14; 
omega_Earth = 7.2921151467e-5;
c = 299792458; 

% Pulling data off from the satellite broadcast and receiver observation
% files
ionospheric_data = read_nav_ion('obs/brdc2970.03n');
ephemeris_data = read_nav_ephemeris('obs/brdc2970.03n');
observation_data = read_obs('obs/2003_297.rnx');

% Setting index numbers for reqiured quantities
GPS_wk_column = 1;
GPS_TOW_column = 2;
PRN_column = 3;
L1_column = 4;
L2_column = 5;
C1_column = 6;

rows = find(observation_data(:,GPS_TOW_column)==min(observation_data(:,GPS_TOW_column)));
PRN_list = observation_data(rows,PRN_column);
GPS_Secs = observation_data(rows,GPS_TOW_column);
GPS_Weeks = observation_data(rows,GPS_wk_column);

% Finding the ephemeris file that is created at a time which is closest to
% our time of observation
[epoch_nav_data,rows] = closest_ephemeris(PRN_list, GPS_Weeks(1),GPS_Secs(1),ephemeris_data)

% Initial observation of receiver position
rec_xyz = obs_initial_pos('obs/2003_297.rnx');
x_rec = rec_xyz(1);
y_rec = rec_xyz(2);
z_rec = rec_xyz(3);
Cb = 0;

% Initialising time and correction quantities
Tt = zeros(length(rows),1);
sat_clk_t_corr = Tt;
satcorr = Tt;
rel_corr = Tt;
C1 = Tt;

%Initialising the required position and range values for satellite and the
%receiver
range = zeros(length(rows),1);
drange = zeros(length(rows),1);
rSat = zeros(length(rows),3);

% Loop for calculating the range with relativity and clock corrections
for ii = 1:length(rows)
    GPS_SOW = epoch_nav_data(ii,17);
    GPS_Week = GPS_Weeks(1);
    Seconds = GPS_Secs(1);
    [R(ii), rSat(ii,1),rSat(ii,2),rSat(ii,3), rel_dt] = sat_geo_range(rec_xyz',GPS_Week,GPS_Secs(1),PRN_list(ii),epoch_nav_data(ii,:),Seconds)
    rel_corr(ii) = rel_dt*c;
    sat_clk_t_corr(ii) = sat_clock_correction(GPS_Week, GPS_Secs(1), PRN_list(ii), epoch_nav_data(ii,:));
    satcorr(ii) = sat_clk_t_corr(ii)*c;
    C1(ii) = observation_data(ii,C1_column);
    % Difference between computed and observed ranges which is to be used
    % in LS and other weighted LS (i.e. EKF, UKF) solutions
    drange(ii) = C1(ii)-R(ii)+satcorr(ii)+rel_corr(ii);
end

% Initialising the system matrix for LS solution   
A = zeros(length(rows),4);

while abs(dx)>1e-5 && abs(dy)>1e-5 && abs(dz)>1e-5 % precision limit for the solution
    
    for ii = 1:length(rows)
        A(ii,:) = [(-(rSat(ii,1)-x_rec)/drange(ii)) (-(rSat(ii,2)-y_rec)/drange(ii)) (-(rSat(ii,3)-z_rec)/drange(ii)) 1];
    end
    
    % LS solution algorithm
    delta = (A'*A)\A'*drange;
    Cb = Cb + delta(4);
    dx = delta(1);
    x_rec = x_rec + dx;
    dy = delta(2);
    y_rec = y_rec +dy;
    dz = delta(3);
    z_rec = z_rec+dz;

    pos_rec = [x_rec,y_rec,z_rec];
end

% Final position of the receiver
fprintf('receiver position: \n')
receiver_position = pos_rec