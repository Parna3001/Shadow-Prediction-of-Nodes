times = cell(1, 10); % Preallocate array for time values
    energies = zeros(1, 10); % Preallocate array for energy values

    fprintf('Date: ');
    date = input('', 's');
    fprintf('Object Height: ');
    height = input('', 's');
    h = str2double(height);
    fprintf('Area of one node: ');
    area = input('', 's');
    ar = str2double(area);
    fprintf('Length: ');
    length = input('', 's');
    l = str2double(length);
    fprintf('Number of nodes: ');
    node = input('', 's');
    n = str2double(node);

for i = 1:10
    fprintf('Enter details for iteration %d:\n', i);
    fprintf('Enter time: ');
    time_input = input('', 's');
    times{i} = time_input;
    energy = enhancedProEnergy(date, time_input, ar, l, n,h);
    energies(i) = energy;
end

% Plot the graph
plot(1:10, energies, 'bo-');
set(gca, 'XTickLabel', times); % Set time labels on x-axis
xlabel('Time');
ylabel('Harvested Energy');
title('Harvested Energy vs. Time');
    


function e = enhancedProEnergy(date,time,ar, l, n,h)

    ha = hourAngle(time);
    disp(['Hour Angle: ', num2str(ha)]);

    d = declination(date, time);
    disp(['Declination Angle: ', num2str(d)]);

    e = elevation(date, time);
    disp(['Elevation Angle: ', num2str(e)]);

    a = azimuth(date, time);
    disp(['Azimuth Angle: ', num2str(a)]);

    z = Zenith(date, time);
    disp(['Zenith Angle: ', num2str(z)]);

    j = jValue(date);
    disp(['Nth day of the year: ', num2str(j)]);

    s = shadowLength(e, h);
    disp(['Shadow Length: ', num2str(s)]);

    dn = directNormal(date,time);
    disp(['Direct Normal Irradiance: ', num2str(dn)]);

    dns = DNIs(date, ar, l, n,time);
    disp(['Direct Normal Irradiance for shaded area: ', num2str(dns)]);

    dhs = calculateDHI(date,time);
    % disp(['Diffuse Horizontal Irradiance for shaded area: ', num2str(dhs)]);
   energy = dhs ;
   disp(['Harvested Energy Short term: ', num2str(energy)]);
   e = energy;
end

% Hour Angle
function h = hourAngle(time)
    t = datetime(time, 'InputFormat', 'HH:mm');
    ho = hour(t);
    h = 15 * (ho - 12);
end

% Day count from 1st January
function count = countDays(date)
    year = str2double(date(1:4));
    month = str2double(date(6:7));
    day = str2double(date(9:end));
    days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    if mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0)
        days(2) = 29;
    end
    count = sum(days(1:month-1)) + day;
end

% Declination Angle
function d = declination(date,time)
    t = datetime(time, 'InputFormat', 'HH:mm');
    ho = hour(t);
    p = ho / 24;
    N = countDays(date) + p; % Extract the day of the year directly using day()
    angle_deg = single((360 * (N + 10)) / 365); % Convert to single precision
    a1 = (cos(deg2rad(angle_deg))) * (-23.44);
    d = a1;
end

% Elevation Angle
function e = elevation(date, time)
    Elevation = asin(((sin(deg2rad(declination(date, time)))) * (sin(deg2rad(22.6)))) + ((cos(deg2rad(22.6))) * (cos(deg2rad(declination(date,time)))) * cos(deg2rad(hourAngle(time)))));
    e = rad2deg(Elevation);
end

% Azimuth Angle
function A = azimuth(date, time)
    Azimuth = acos(((sin(deg2rad(declination(date, time)))) * (cos(deg2rad(22.6)))) - ...
        ((sin(deg2rad(22.6))) * (cos(deg2rad(declination(date,time)))) * cos(deg2rad(hourAngle(time)))) / (cos(deg2rad(elevation(date,time)))));
    if hourAngle(time) > 0
        azimuth = deg2rad(360) - Azimuth;
    else
        azimuth = Azimuth;
    end
    A = rad2deg(azimuth);
end

% Zenith Angle
function z = Zenith(date, time)
    z = (90 - elevation(date, time));
end

% Nth Day of the year
function j = jValue(date)
    d = datetime(date, 'InputFormat', 'yyyy-MM-dd');
    j = day(d, 'dayofyear');
end

% Shadow Length
function s = shadowLength(E, h)
    a = rad2deg(tand(deg2rad(E)));
    s = (h / a);
end
% Direct Normal Irradiance
function DNI = directNormal(date, time)
    % Constants
    n = 0.8;  % Example value for clearness index
    
    % Calculate solar altitude angle
    Elevation = elevation(date, time);
    
    % Calculate solar zenith angle
    ZenithAngle = 90 - Elevation;
    
    % Calculate the clearness index correction factor
    if n < 0.95
        tau_c = (n + 0.05);
    elseif n == 0.95
        tau_c = 1.00;
    end
    
    % Calculate the relative optical air mass
    m = 1 / cos(deg2rad(ZenithAngle));
    
    % Calculate the atmospheric transmittance for direct solar radiation
    tau_atmdir = (exp(-0.65 * m) + exp(-0.95 * m)) / 2.0;
    
    % Calculate the cloud transmittance for direct solar radiation
    if m ~= tau_c
        tau_clouddir = (1 - tau_c) / (m - tau_c);
    end
    j = jValue(date);
    DNI = ((1366.1 * (1 + 0.033 * (rad2deg(cos(deg2rad(((360 * j) / 365)))))))) * tau_atmdir * tau_clouddir;
end
% Direct Normal Irradiance For Shaded Area
function DNs = DNIs(date, s, l, node,time)
    DNs = directNormal(date,time) * (1 - ((s * l) / (l * l * node)));
end

% Diffuse Horizontal Irradiance For Shaded area
function DHs = calculateDHI(date,time)
% Perez model
    A = 0.1787;
    B = 0.5065;
    C = -0.5593;
    D = 0.8760;
    am = max(0, pi/2 - deg2rad(elevation(date,time)));
    ai = asin(sin(am)) * sin(deg2rad(azimuth(date,time)));
    diffuseFactor = A + B * exp(C/am) + D * cos(ai)^2;
    diffuseFactor = diffuseFactor * (1 + cos(deg2rad(9.14))) / 2;
    diffuseFactor = diffuseFactor / directNormal(date,time);
    d = diffuseFactor * (directNormal(date,time)) / 5.06;
    DHI = 5.06 - (directNormal(date,time) * (rad2deg(cos(deg2rad(azimuth(date, time))))));
    DHs = DHI * (1 - d);
end

% 0	33.65	120.85	236.85	   345.27	428.38	411.4	375.08	292.73	183.3	80.15	10.7