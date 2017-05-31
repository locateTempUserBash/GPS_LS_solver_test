function [ionospheric_params]=read_nav_ion(nav_fname)
% This function reads ionospheric parameters from navigation data files
% ending with ".03n"

FID = fopen(nav_fname);

headerend   = [];
headeralpha = [];  alpha = [];
headerbeta  = [];  beta  = [];

while (isempty(headerend) == 1)
   nav_line     = fgetl(FID); 
   headerend = findstr(nav_line,'END OF HEADER');
   
   headeralpha = findstr(nav_line,'ION ALPHA');
   if (isempty(headeralpha) == 0)
       
      [A0, remain] = strtok(nav_line);
      [A1, remain] = strtok(remain);
      [A2, remain] = strtok(remain);
      [A3]         = strtok(remain);
      
      alpha = [str2num(A0) str2num(A1) str2num(A2) str2num(A3)];
      
   end
   
   headerbeta = findstr(nav_line,'ION BETA');
   if (isempty(headerbeta) == 0)
       
      [B0, remain] = strtok(nav_line);
      [B1, remain] = strtok(remain);
      [B2, remain] = strtok(remain);
      [B3]         = strtok(remain);
      
      beta = [str2num(B0) str2num(B1) str2num(B2) str2num(B3)];
      
   end
   
end

ionospheric_params = [alpha beta];

fclose(FID);
end