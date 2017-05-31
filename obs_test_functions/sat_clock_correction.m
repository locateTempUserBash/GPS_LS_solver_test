function [tcorr] = sat_clock_correction(GPS_Weeks, GPS_SOW, PRN, nav_data) 


Af0_col = 21;   %Af0, satellite clock bias (sec)
Af1_col = 22;   %Af1, satellite clock drift (sec/sec)
Af2_col = 23;   %Af2, satellite clock drift rate (sec/sec/sec)
SOW_col = 17;   %Toe, reference time ephemeris (seconds into GPS week)

Af0 = nav_data(Af0_col); 
Af1 = nav_data(Af1_col); 
Af2 = nav_data(Af2_col);

t_eph = nav_data(SOW_col);
dt = GPS_SOW - t_eph;

tcorr = Af0 + Af1*(dt) + Af2*(dt)^2;

end 



