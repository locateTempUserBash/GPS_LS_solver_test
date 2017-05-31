function [obs,Gregorian_time] = read_obs(obs_fname)

% This function reads L1, L2 and C1 data from rinex observation files
% with format "rnx"

obs = [];
total_satellites=[];
obs_data=[];
Gregorian_time=zeros(1,6);
FID = fopen(obs_fname);

% The first 13 lines are the Header

for i = 1:14
    obs_line = fgetl(FID);
   
end

% Drawing time and observation data

k=0;j=1;
while 1
    
    obs_line = obs_line(2:length(obs_line));
    F = sscanf(obs_line,'%f');
    
    if (size(F)==0), break, end
    
    % Load GPS Gregorian time into variables
    
        Y = F(1) + 2000;
        M = F(2);
        D = F(3);
        H = F(4);
        min = F(5);
        sec = F(6);
        Gregorian_time(j,:) = [Y M D H min sec];
        utc_time = [Y M D H min sec];
        
        % Convert GPS Gregorian time to GPS week and GPS TOW
        [GPS_wk, GPS_TOW] = date2GPSTime(utc_time) ;
        
        % Number of Satellites and the PRN list
        total_satellites = [total_satellites;F(8)];
        PRN = F(9:length(F));
        
        % Next block of lines have the observation data
        obs_line = fgetl(FID);
        for n=1:length(PRN)
            % Line is scaned for each satellite
            obs_line = obs_line(3:length(obs_line));
            F = sscanf(obs_line,'%f');
            L1 = F(1);
            L2 = F(2);
            C1 = F(3);
            PRN_no = PRN(n);
            obs(k+n,:) = [GPS_wk GPS_TOW PRN_no L1 L2 C1];
            
            obs_line = fgetl(FID); % next line is irrelevant data, but
            obs_line = fgetl(FID); % the following line has L1, L2 and C1
                                   % data
        end
        
        k = k+n-1;
        j = j+1;
        
       

        
end

fclose(FID);

end

            
       