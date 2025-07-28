% Space-Time Spectral Analysis for CESM2 CAM data at 300 hPa
% Historical and RCP8.5: Latitude–phase speed spectrum of eddy momentum flux convergence

close all; clear; clc;

RADIUS = 6.371e6;
DAY = 86400;

% Set path to data files
histfile = 'b40.20th.track1.1deg.cam2.h0.1961-1980.1x1.nc';
rcpfile  = 'b40.rcp8_5.1deg.cam2.h0.2081-2100.1x1.nc';

% --- Choose pressure level index (lev_idx) for 300 hPa
lev = ncread(histfile,'lev');
[~,lev_idx] = min(abs(lev-300)); % Correct index for 300 hPa

do_compute = 1; % Set to 1 to compute, 0 to only plot/load

if do_compute > 0
    %% READ DATA
    % --- Historical
    lat = ncread(histfile, 'lat');    % 180
    lon = ncread(histfile, 'lon');    % 360
    U_hist = ncread(histfile, 'U');   % (lon, lat, lev, time)
    V_hist = ncread(histfile, 'V');
    % select one level
    U_hist = squeeze(U_hist(:,:,lev_idx,:)); % (lon,lat,time)
    V_hist = squeeze(V_hist(:,:,lev_idx,:));
    U_hist = permute(U_hist, [3,2,1]); % (time,lat,lon)
    V_hist = permute(V_hist, [3,2,1]);

    % --- RCP85
    U_rcp  = ncread(rcpfile, 'U');
    V_rcp  = ncread(rcpfile, 'V');
    U_rcp  = squeeze(U_rcp(:,:,lev_idx,:));
    V_rcp  = squeeze(V_rcp(:,:,lev_idx,:));
    U_rcp  = permute(U_rcp, [3,2,1]); % (time,lat,lon)
    V_rcp  = permute(V_rcp, [3,2,1]);

    rad = lat/180*pi;

    % --- Partition data into years (20 years, 12 months each)
    nyr = 20; nmon = 12; nlat = length(lat); nlon = length(lon);
    uu_yy_hist = zeros(nyr, nmon, nlat, nlon);
    vv_yy_hist = zeros(nyr, nmon, nlat, nlon);
    uu_yy_rcp  = zeros(nyr, nmon, nlat, nlon);
    vv_yy_rcp  = zeros(nyr, nmon, nlat, nlon);
    for yy = 1:nyr
        t_idx = (yy-1)*nmon + (1:nmon);
        uu_yy_hist(yy,:,:,:) = U_hist(t_idx,:,:);
        vv_yy_hist(yy,:,:,:) = V_hist(t_idx,:,:);
        uu_yy_rcp(yy,:,:,:)  = U_rcp(t_idx,:,:);
        vv_yy_rcp(yy,:,:,:)  = V_rcp(t_idx,:,:);
    end

    hann_win = hann(nmon);

    % --- Space-time analysis (historical)
    for yy=1:nyr
        u = zeros(nlon, nmon, nlat);
        v = zeros(nlon, nmon, nlat);
        for i=1:nlon
            for j=1:nlat
                for t=1:nmon
                    u(i,t,j)=squeeze(uu_yy_hist(yy,t,j,i))*hann_win(t);
                    v(i,t,j)=squeeze(vv_yy_hist(yy,t,j,i))*hann_win(t);
                end
            end
        end
        [East, West]=space_time_new(u,v); % (k, freq, lat)
        [nk, nf, nlat_eff] = size(East);
        % Preallocate on first loop
        if yy==1
            K_e_m_yy_hist = zeros(nyr, nk, nf, nlat_eff);
            K_w_m_yy_hist = zeros(nyr, nk, nf, nlat_eff);
        end
        % Multiply for cos(lat)^2 and compute meridional derivative
        K_e_mm_cos_hist = zeros(nk, nf, nlat_eff);
        K_w_mm_cos_hist = zeros(nk, nf, nlat_eff);
        for ff=1:nf
            for jj=1:nlat_eff
                K_e_mm_cos_hist(:,ff,jj)=East(:,ff,jj)*(cos(rad(jj))^2);
                K_w_mm_cos_hist(:,ff,jj)=West(:,ff,jj)*(cos(rad(jj))^2);
            end
        end
        dK_e_m=dy_dx(K_e_mm_cos_hist,rad,3)/RADIUS*DAY;
        dK_w_m=dy_dx(K_w_mm_cos_hist,rad,3)/RADIUS*DAY;
        for ff=1:nf
            for jj=1:nlat_eff
                K_e_m_yy_hist(yy,:,ff,jj)=-dK_e_m(:,ff,jj)./(cos(rad(jj))^2);
                K_w_m_yy_hist(yy,:,ff,jj)=-dK_w_m(:,ff,jj)./(cos(rad(jj))^2);
            end
        end
    end
    % --- Space-time analysis (RCP)
    for yy=1:nyr
        u = zeros(nlon, nmon, nlat);
        v = zeros(nlon, nmon, nlat);
        for i=1:nlon
            for j=1:nlat
                for t=1:nmon
                    u(i,t,j)=squeeze(uu_yy_rcp(yy,t,j,i))*hann_win(t);
                    v(i,t,j)=squeeze(vv_yy_rcp(yy,t,j,i))*hann_win(t);
                end
            end
        end
        [East, West]=space_time_new(u,v); % (k, freq, lat)
        [nk, nf, nlat_eff] = size(East);
        if yy==1
            K_e_m_yy_rcp = zeros(nyr, nk, nf, nlat_eff);
            K_w_m_yy_rcp = zeros(nyr, nk, nf, nlat_eff);
        end
        K_e_mm_cos_rcp = zeros(nk, nf, nlat_eff);
        K_w_mm_cos_rcp = zeros(nk, nf, nlat_eff);
        for ff=1:nf
            for jj=1:nlat_eff
                K_e_mm_cos_rcp(:,ff,jj)=East(:,ff,jj)*(cos(rad(jj))^2);
                K_w_mm_cos_rcp(:,ff,jj)=West(:,ff,jj)*(cos(rad(jj))^2);
            end
        end
        dK_e_m=dy_dx(K_e_mm_cos_rcp,rad,3)/RADIUS*DAY;
        dK_w_m=dy_dx(K_w_mm_cos_rcp,rad,3)/RADIUS*DAY;
        for ff=1:nf
            for jj=1:nlat_eff
                K_e_m_yy_rcp(yy,:,ff,jj)=-dK_e_m(:,ff,jj)./(cos(rad(jj))^2);
                K_w_m_yy_rcp(yy,:,ff,jj)=-dK_w_m(:,ff,jj)./(cos(rad(jj))^2);
            end
        end
    end

    %=== smooth and zonal mean
    K_e_sm_yy_hist = zeros(size(K_e_m_yy_hist));
    K_w_sm_yy_hist = zeros(size(K_w_m_yy_hist));
    K_e_sm_yy_rcp  = zeros(size(K_e_m_yy_rcp));
    K_w_sm_yy_rcp  = zeros(size(K_w_m_yy_rcp));
    K_e_m_yy_hist(:,1,:,:) = 0;   K_w_m_yy_hist(:,1,:,:) = 0;
    K_e_m_yy_rcp(:,1,:,:)  = 0;   K_w_m_yy_rcp(:,1,:,:)  = 0;
    for nn=1:nyr
        for j=1:nlat_eff
            % Concatenate over wavenumbers 2:end only, so size matches
            A_hist = squeeze(K_w_m_yy_hist(nn,2:end,:,j));   % (nk-1, nf)
            B_hist = squeeze(K_e_m_yy_hist(nn,2:end,:,j));   % (nk-1, nf)
            
            % Check if we have valid data before smoothing
            if ~isempty(A_hist) && ~isempty(B_hist) && size(A_hist,2) > 0 && size(B_hist,2) > 0
                KK_hist = smooth([fliplr(A_hist), B_hist],3,'g'); % (nk-1, 2*nf)
                % Extract the eastward component (second half)
                nf_cols = size(B_hist,2);
                if size(KK_hist,2) >= nf_cols
                    K_e_sm_yy_hist(nn,2:end,:,j) = KK_hist(:,nf_cols+1:end);  % eastward
                    K_w_sm_yy_hist(nn,2:end,:,j) = fliplr(KK_hist(:,1:nf_cols)); % westward
                else
                    % Fallback: use original data if smoothing fails
                    K_e_sm_yy_hist(nn,2:end,:,j) = B_hist;
                    K_w_sm_yy_hist(nn,2:end,:,j) = A_hist;
                end
            else
                % Set to zero if no valid data
                K_e_sm_yy_hist(nn,2:end,:,j) = 0;
                K_w_sm_yy_hist(nn,2:end,:,j) = 0;
            end
            
            A_rcp = squeeze(K_w_m_yy_rcp(nn,2:end,:,j));
            B_rcp = squeeze(K_e_m_yy_rcp(nn,2:end,:,j));
            
            % Same fix for RCP data
            if ~isempty(A_rcp) && ~isempty(B_rcp) && size(A_rcp,2) > 0 && size(B_rcp,2) > 0
                KK_rcp = smooth([fliplr(A_rcp), B_rcp],3,'g');
                nf_cols = size(B_rcp,2);
                if size(KK_rcp,2) >= nf_cols
                    K_e_sm_yy_rcp(nn,2:end,:,j) = KK_rcp(:,nf_cols+1:end);
                    K_w_sm_yy_rcp(nn,2:end,:,j) = fliplr(KK_rcp(:,1:nf_cols));
                else
                    K_e_sm_yy_rcp(nn,2:end,:,j) = B_rcp;
                    K_w_sm_yy_rcp(nn,2:end,:,j) = A_rcp;
                end
            else
                K_e_sm_yy_rcp(nn,2:end,:,j) = 0;
                K_w_sm_yy_rcp(nn,2:end,:,j) = 0;
            end
        end
    end

    uzm_hist = squeeze(mean(mean(uu_yy_hist,4),2)); % (years, lat)
    uzm_rcp  = squeeze(mean(mean(uu_yy_rcp,4),2));

    save('uv_DJFM_CESM2_hist.mat','K_e_sm_yy_hist','K_w_sm_yy_hist','uzm_hist','lat');
    save('uv_DJFM_CESM2_rcp85.mat','K_e_sm_yy_rcp','K_w_sm_yy_rcp','uzm_rcp','lat');
end

%% Plotting and further processing
load('uv_DJFM_CESM2_hist.mat','K_e_sm_yy_hist','K_w_sm_yy_hist','uzm_hist','lat');
uc_hist = squeeze(mean(uzm_hist,1))';
K_e_m = squeeze(mean(K_e_sm_yy_hist,1));
K_w_m = squeeze(mean(K_w_sm_yy_hist,1));

[nk, nf, nlat_eff] = size(K_w_m);
wavenum = 0:(nk-1);
freq = 0.5*(0:(nf-1))/(nf-1);
mm = 2:nk;
wavenum = wavenum(mm);
K_e = K_e_m(mm,:,:);
K_w = K_w_m(mm,:,:);
nk2 = length(wavenum);  nj = nlat_eff;
kk = 1:nk2; jj=1:nj;
c_max=50;  c_min = -10;   ca=0.0:0.1:c_max;   cc=1:length(ca);

for pos=1:nj
    [K_e_new, K_w_new] = convert_fk_to_ck(freq, wavenum(kk), cos(lat(pos)/180.*pi), ca, squeeze(K_e(:, :, pos)), squeeze(K_w(:, :, pos)));
    K_e_lat_hist(:,pos)=sum(K_e_new,2);
    K_w_lat_hist(:,pos)=sum(K_w_new,2);
end

clear K_e_sm_yy_hist K_w_sm_yy_hist uzm_hist

load('uv_DJFM_CESM2_rcp85.mat','K_e_sm_yy_rcp','K_w_sm_yy_rcp','uzm_rcp','lat');
uc_rcp85 = squeeze(mean(uzm_rcp,1))';
K_e_m = squeeze(mean(K_e_sm_yy_rcp,1));
K_w_m = squeeze(mean(K_w_sm_yy_rcp,1));
K_e = K_e_m(mm,:,:);
K_w = K_w_m(mm,:,:);

for pos=1:nj
    [K_e_new, K_w_new] = convert_fk_to_ck(freq, wavenum(kk), cos(lat(pos)/180.*pi), ca, squeeze(K_e(:, :, pos)), squeeze(K_w(:, :, pos)));
    K_e_lat_rcp85(:,pos)=sum(K_e_new,2);
    K_w_lat_rcp85(:,pos)=sum(K_w_new,2);
end

save('co-spectra_cesm2.mat','K_w_lat_hist','K_w_lat_rcp85','K_e_lat_hist','K_e_lat_rcp85','uc_hist','uc_rcp85');

figure('pos',[10 10 700 700]);
step = 0.01;
colorb = -5*step:step:step*5;
[cc1,hh1]=contourf(-ca(cc),lat(jj),1.25*squeeze(K_w_lat_rcp85(cc,jj)-K_w_lat_hist(cc,jj))',colorb); hold on;
[cc2,hh2]=contourf(ca(cc),lat(jj),1.25*squeeze(K_e_lat_rcp85(cc,jj)-K_e_lat_hist(cc,jj))',colorb);
caxis([min(colorb),max(colorb)]);
colormap(jet); set([hh1,hh2],'edgecolor','none');
step = 0.02;
[c1,h1]=contour(-ca(cc),lat(jj),squeeze(K_w_lat_hist(cc,jj))',step:step:step*12,'k-','linewidth',2); hold on;
[c2,h2]=contour(ca(cc), lat(jj),squeeze(K_e_lat_hist(cc,jj))',step:step:step*12,'k-','linewidth',2);
axis([c_min,c_max,-80,80]);
set([h1,h2],'linewidth',2.0)
xlabel('angular phase speed (m/s)','fontsize',20);
ylabel('latitude(deg)','fontsize',20);
set(gca,'xtick',c_min:5:c_max,'ytick',-80:20:80,'ylim',[-80,80],'fontsize',20);
line([0,0],[-80,80],'linestyle','-.','color','k');
plot(uc_rcp85./(cos(lat*pi/180)),lat,'r','linewidth',2);hold on;
plot(uc_hist./(cos(lat*pi/180)),lat,'k','linewidth',2);

% --- Helper functions must be in separate .m files in your MATLAB path:
% - space_time_new.m
% - smooth.m
% - dy_dx.m
% - convert_fk_to_ck.m