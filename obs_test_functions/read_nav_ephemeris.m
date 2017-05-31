function [gps_ephemeris] = read_nav_ephemeris(nav_fname)
% This function reads gps ephemeris data from navigation files ending with
% ".03n" suffix

FID = fopen(nav_fname);
headerend = [];

while (isempty(headerend) == 1)
   nav_line     = fgetl(FID); 
   headerend = findstr(nav_line,'END OF HEADER');
end
j = 1;


while 1
    
    nav_line = fgetl(FID);
    
    % If the next line is not a character then the end of the file has been
    %   reached and the while loop is exited
    if ~ischar(nav_line), break, end

        % Line 1
        prn = str2num(nav_line(1:2));   % PRN number of satellite    
        Af0 = str2num(nav_line(23:41)); % clock bias  (s)
        Af1 = str2num(nav_line(42:60)); % clock drift (s/s)
        Af2 = str2num(nav_line(61:79)); % clock drift rate (s/s/s)
        %-----
        
        % Line 2
        nav_line = fgetl(FID); % Read in second line of satellite ephemeris
        IODE = str2num(nav_line(4:22)) ;   % Issue of Data (Ephemeris)       
        Crs  = str2num(nav_line(23:41));   % Amplitude of the Sine Harmonic Correction 
                                           %  Term to the orbit radius
      
        delta_n= str2num(nav_line(42:60)); % Mean Motion Difference from Computed Value, rad/s
        
        M0   = str2num(nav_line(61:79));   % Mean Anomaly at Reference Time, rad
        %-----
            
        % Line 3
        nav_line = fgetl(FID); % Read in third line of satellite ephemeris
        Cuc = str2num(nav_line(4:22)) ;    % Amplitude of the Cosine Harmonic Correction
                                        %  Term to the Argument of Latitude
        ecc = str2num(nav_line(23:41)) ;   % Eccentricity
        
        Cus = str2num(nav_line(42:60));    % Amplitude of the Sine Harmonic Correction
                                        %  Term to the Argument of Latitude
       
        sqrt_a = str2num(nav_line(61:79)); % Square root of the semi-major Axis, m^0.5
        %-----
        
        % Line 4
        nav_line = fgetl(FID); % Read in fourth line of satellite ephemeris
        Toe = str2num(nav_line(4:22));     % Reference Time Ephemeris (sec into GPS week)
        Cic = str2num(nav_line(23:41));    % Amplitude of the Cosine Harmonic Correction
                                        %  Term to the Angle of Inclination
        Loa = str2num(nav_line(42:60));    % Longitude of Ascending Node of Orbit Plane
                                        %  at Weekly Epoch, rad
      
        Cis = str2num(nav_line(61:79)) ;   % Amplitude of the Sine Harmonic Correction
                                        %  Term to the Angle of Inclination
        %-----
        
        % Line 5 
        nav_line = fgetl(FID); % Read in fifth line of satellite ephemeris
        
        incl = str2num(nav_line(4:22));    % Inclination Angle at Reference Time, rad
        
        Crc  = str2num(nav_line(23:41));   % Amplitude of the Cosine Harmonic Correction
                                        %  Term to the Orbit Radius
     
        perigee = str2num(nav_line(42:60));% Argument of Perigee, rad
        
        ra_rate = str2num(nav_line(61:79));% Rate of Change of Right Ascension, rad/s
       
        %Line 6
        nav_line = fgetl(FID); % Read in sixth line of satellite ephemeris
        
        i_rate = str2num(nav_line(4:22));   % Rate of change of inclination angle, rad/s
        
        %str = tline(23:41);             % codes on L2 channel (unecessary)
       
        GPS_week = str2num(nav_line(42:60));% GPS Week Number (to go with Toe)
        
        %str   = tline(61:79);           % L2 flag
        
        %-----
        
        % Line 7
        nav_line = fgetl(FID); % Includes: SV accuracy, SV health, TGD, IODC
        
        health = str2num(nav_line(23:41)); % Satellite health (0.00 = usable)
        %-----
        
        
        % Line 8
        nav_line = fgetl(FID); % Read in eighth line of satellite ephemeris
        
        Toc = Toe; % Time of clock
        
        
        
        
        ephemeris_vector = [prn M0 delta_n ecc...
            sqrt_a Loa incl perigee ...
            ra_rate i_rate Cuc Cus Crc Crs Cic Cis Toe IODE ...
            GPS_week Toc Af0 Af1 Af2 0 health];
        vector_size = length(ephemeris_vector);
        gps_ephemeris(j,1:vector_size)=ephemeris_vector;
        
        j = j + 1;   
        
end


fclose(FID);
end
