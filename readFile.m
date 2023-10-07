function [opened, vehID, GST, date_time, R, V] = fileRead(fileName)
    
    opened = 1;
    vehID = NaN;
    R = NaN;
    V = NaN;
    % Open the file for reading
    try
        fileID = fopen(fileName, 'r');

        year = fscanf(fileID, '$ YEAR %d', 1); % Reading in year variable (line 1)
        month_num = fscanf(fileID, '   MNTH %d', 1); % Reading in month variable (line 1)
        day = fscanf(fileID, '  DAY %d\n', 1); % Reading in day variable (line 1)

        % If statement to give abbreviated month 
        if month_num == 1
            month = 'JAN';
        elseif month_num == 2
            month = 'FEB';
        elseif month_num == 3
            month = 'MAR';
        elseif month_num == 4
            month = 'APR';
        elseif month_num == 5
            month = 'MAY';
        elseif month_num == 6
            month = 'JUN';
        elseif month_num == 7
            month = 'JUL';
        elseif month_num == 8
            month = 'AUG';
        elseif month_num == 9
            month = 'SEP';
        elseif month_num == 10
            month = 'OCT';
        elseif month_num == 11
            month = 'NOV';
        else 
            month = 'DEC';
        end

        date = [month,"/",day,"/",year]; % Concatenating variables into date format

        hour = fscanf(fileID, '$ HR %d', 1); % Reading in hour variable (line 2)
        minute = fscanf(fileID, '   MIN %d', 1); % Reading in minute variable (line 2)
        seconds = fscanf(fileID, '  SEC %d', 1); % Reading in seconds variable (line 2)

        GST = ((hour*3600)+(minute*60)+seconds); % Converting time into seconds

        time = [hour,":",minute,":",seconds]; % Concatenating time variables

        date_time = [date," ",time]; % Concatenating date and time arrays

        vehID = fscanf(fileID, '$ VEHID %f\n', 1); % line 4
    
        % skip fifth line
        fgetl(fileID);

        Rx = fscanf(fileID, '$ %f\n', 1)/1000; % (km) line 6
        Ry = fscanf(fileID, '$ %f\n', 1)/1000; % (km) line 7
        Rz = fscanf(fileID, '$ %f\n', 1)/1000; % (km) line 8
        R = [Rx, Ry, Rz];

        Vx = fscanf(fileID, '$ %f\n', 1)/1000; % (km/s) line 9
        Vy = fscanf(fileID, '$ %f\n', 1)/1000; % (km/s) line 10
        Vz = fscanf(fileID, '$ %f\n', 1)/1000; % (km/s) line 11
        V = [Vx, Vy, Vz];

        % Close the file
        fclose(fileID);
    catch exception
        opened = 0;
    end

end