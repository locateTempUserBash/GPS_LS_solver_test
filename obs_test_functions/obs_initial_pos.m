function [xyz] = obs_initial_pos(obs_fname)
headerend = [];
headerxyz = [];
FID = fopen(obs_fname);

while (isempty(headerend) == 1)
    obs_line = fgetl(FID);
    headerend = findstr(obs_line,'END OF HEADER');
    headerxyz = findstr(obs_line,'APPROX POSITION XYZ');
    if (isempty(headerxyz) == 0)
        F = sscanf(obs_line,'%f');
        x = F(1); y = F(2); z = F(3);
        xyz = [x; y; z];
    end
end

fclose(FID);
end