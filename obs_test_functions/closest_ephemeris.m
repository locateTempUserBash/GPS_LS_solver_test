function [ephemData,rownums] = closest_ephemeris(PRN_No, GPS_Weeks, GPS_SOW, gps_ephemeris)

rownums = find(gps_ephemeris(:,17)<=GPS_SOW & ismember(gps_ephemeris(:,1),PRN_No) & gps_ephemeris(:,19)==GPS_Weeks);
ephemData = gps_ephemeris(rownums,:);
