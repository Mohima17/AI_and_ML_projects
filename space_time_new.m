function [East, West] = space_time_new(u, v)
% Space-time spectral analysis function
% Input: u, v are (nlon, ntime, nlat) arrays
% Output: East, West are (wavenumber, frequency, latitude) arrays

[nlon, ntime, nlat] = size(u);

% Initialize output arrays
nk = nlon;
nf = ntime;
East = zeros(nk, nf, nlat);
West = zeros(nk, nf, nlat);

% Perform space-time FFT for each latitude
for j = 1:nlat
    % Extract u and v at this latitude
    u_lat = squeeze(u(:, :, j));  % (nlon, ntime)
    v_lat = squeeze(v(:, :, j));  % (nlon, ntime)
    
    % 2D FFT in space and time
    u_fft = fft2(u_lat);
    v_fft = fft2(v_lat);
    
    % Calculate power spectral density
    u_power = abs(u_fft).^2;
    v_power = abs(v_fft).^2;
    
    % Separate eastward and westward components
    % This is a simplified approach - you may need to adjust based on your specific needs
    for k = 1:nk
        for f = 1:nf
            % Simple separation based on phase relationship
            % Eastward: positive frequency with positive wavenumber
            % Westward: positive frequency with negative wavenumber
            if k <= nk/2
                East(k, f, j) = u_power(k, f) + v_power(k, f);
                West(k, f, j) = 0;
            else
                East(k, f, j) = 0;
                West(nk - k + 1, f, j) = u_power(k, f) + v_power(k, f);
            end
        end
    end
end

return