close all; clear; 
%param; 
cd /workspace/
 
RADIUS=6.371e6;
DAY=86400;
 
do_compute=1;
 
if do_compute>0 
%=========================== read data ================================
% Read historical data (1961-1980)
hist_file = 'b40.20th.track1.1deg.cam2.h0.1961-1980.1x1.nc';
% Read future data (2081-2100) 
rcp_file = 'b40.rcp8_5.1deg.cam2.h0.2081-2100.1x1.nc';

% Read coordinates
lat = ncread(hist_file, 'lat');
lon = ncread(hist_file, 'lon');
lev = ncread(hist_file, 'lev');

% Find 300 hPa level (index 16 as specified)
lev_300_idx = 16;

% Read U and V wind components at 300 hPa level
% Historical period (1961-1980)
uu_hist = squeeze(ncread(hist_file, 'U', [1 lev_300_idx 1 1], [Inf 1 Inf Inf]));
vv_hist = squeeze(ncread(hist_file, 'V', [1 lev_300_idx 1 1], [Inf 1 Inf Inf]));

% Future period (2081-2100)
uu_rcp = squeeze(ncread(rcp_file, 'U', [1 lev_300_idx 1 1], [Inf 1 Inf Inf]));
vv_rcp = squeeze(ncread(rcp_file, 'V', [1 lev_300_idx 1 1], [Inf 1 Inf Inf]));

% Rearrange dimensions: [time, lat, lon] -> [lon, time, lat]
uu_hist = permute(uu_hist, [3, 1, 2]);
vv_hist = permute(vv_hist, [3, 1, 2]);
uu_rcp = permute(uu_rcp, [3, 1, 2]);
vv_rcp = permute(vv_rcp, [3, 1, 2]);

rad = lat/180.*pi;

% Use all monthly data (240 months = 20 years)
num_years = 20;
months_per_year = 12;

% Initialize arrays for yearly data
uu_yy_hist = zeros(num_years, months_per_year, size(uu_hist,3), size(uu_hist,1));
vv_yy_hist = zeros(num_years, months_per_year, size(vv_hist,3), size(vv_hist,1));
uu_yy_rcp = zeros(num_years, months_per_year, size(uu_rcp,3), size(uu_rcp,1));
vv_yy_rcp = zeros(num_years, months_per_year, size(vv_rcp,3), size(vv_rcp,1));

% Reshape data into yearly chunks (12 months per year)
for yy = 1:num_years
    start_idx = (yy-1)*months_per_year + 1;
    end_idx = yy*months_per_year;
    uu_yy_hist(yy,:,:,:) = uu_hist(:,start_idx:end_idx,:);
    vv_yy_hist(yy,:,:,:) = vv_hist(:,start_idx:end_idx,:);
    uu_yy_rcp(yy,:,:,:) = uu_rcp(:,start_idx:end_idx,:);
    vv_yy_rcp(yy,:,:,:) = vv_rcp(:,start_idx:end_idx,:);
end

%=========================== space-time analysis ================================

% Apply Hanning window (12 months)
hann_win = hann(months_per_year);

% Process historical data
for yy = 1:num_years   % beginning of loop for historical data
 
% ------------------ applying a hanning window ----------------------
for i = 1:size(uu_yy_hist,4)
    for j = 1:size(uu_yy_hist,3)
        for t = 1:size(uu_yy_hist,2)
            u_hist(i,t,j) = uu_yy_hist(yy,t,j,i)*hann_win(t);
            v_hist(i,t,j) = vv_yy_hist(yy,t,j,i)*hann_win(t);
        end
    end
end
 
% ---------------------- space-time spectrum -------------------------
[East_hist, West_hist] = space_time_new(u_hist, v_hist);
 
% --------------------- momentum flux divergence ---------------------
for ff = 1:size(East_hist,2)
    for jj = 1:size(East_hist,3)
        K_e_mm_cos_hist(ff,jj,:) = East_hist(:,ff,jj)*(cos(rad(jj))^2);
        K_w_mm_cos_hist(ff,jj,:) = West_hist(:,ff,jj)*(cos(rad(jj))^2);
    end
end

dK_e_m_hist = dy_dx(K_e_mm_cos_hist,rad,2)/RADIUS*DAY;
dK_w_m_hist = dy_dx(K_w_mm_cos_hist,rad,2)/RADIUS*DAY;

for ff = 1:size(East_hist,2)
    for jj = 1:size(East_hist,3)
        K_e_m_yy_hist(yy,ff,jj,:) = -dK_e_m_hist(ff,jj,:)/(cos(rad(jj))^2);
        K_w_m_yy_hist(yy,ff,jj,:) = -dK_w_m_hist(ff,jj,:)/(cos(rad(jj))^2);
    end
end
 
end           % end of loop for historical data

% Process RCP8.5 data
for yy = 1:num_years   % beginning of loop for RCP8.5 data
 
% ------------------ applying a hanning window ----------------------
for i = 1:size(uu_yy_rcp,4)
    for j = 1:size(uu_yy_rcp,3)
        for t = 1:size(uu_yy_rcp,2)
            u_rcp(i,t,j) = uu_yy_rcp(yy,t,j,i)*hann_win(t);
            v_rcp(i,t,j) = vv_yy_rcp(yy,t,j,i)*hann_win(t);
        end
    end
end
 
% ---------------------- space-time spectrum -------------------------
[East_rcp, West_rcp] = space_time_new(u_rcp, v_rcp);
 
% --------------------- momentum flux divergence ---------------------
for ff = 1:size(East_rcp,2)
    for jj = 1:size(East_rcp,3)
        K_e_mm_cos_rcp(ff,jj,:) = East_rcp(:,ff,jj)*(cos(rad(jj))^2);
        K_w_mm_cos_rcp(ff,jj,:) = West_rcp(:,ff,jj)*(cos(rad(jj))^2);
    end
end

dK_e_m_rcp = dy_dx(K_e_mm_cos_rcp,rad,2)/RADIUS*DAY;
dK_w_m_rcp = dy_dx(K_w_mm_cos_rcp,rad,2)/RADIUS*DAY;

for ff = 1:size(East_rcp,2)
    for jj = 1:size(East_rcp,3)
        K_e_m_yy_rcp(yy,ff,jj,:) = -dK_e_m_rcp(ff,jj,:)/(cos(rad(jj))^2);
        K_w_m_yy_rcp(yy,ff,jj,:) = -dK_w_m_rcp(ff,jj,:)/(cos(rad(jj))^2);
    end
end
 
end           % end of loop for RCP8.5 data

%=========================== smooth ================================
K_e_sm_yy_hist = zeros(size(K_e_m_yy_hist));
K_w_sm_yy_hist = zeros(size(K_w_m_yy_hist));
K_e_sm_yy_rcp = zeros(size(K_e_m_yy_rcp));
K_w_sm_yy_rcp = zeros(size(K_w_m_yy_rcp));

% removing stationary waves
K_e_m_yy_hist(:,1,:,:) = 0;
K_w_m_yy_hist(:,1,:,:) = 0;
K_e_m_yy_rcp(:,1,:,:) = 0;
K_w_m_yy_rcp(:,1,:,:) = 0;

% Smooth historical data
for nn = 1:num_years
    for j = 1:length(lat)
        for i = 1:size(K_e_m_yy_hist,4)
            if size(K_w_m_yy_hist,2) >= 3 && size(K_e_m_yy_hist,2) >= 3
                KK_hist = smooth([fliplr(squeeze(K_w_m_yy_hist(nn,2:end,j,i))'), squeeze(K_e_m_yy_hist(nn,1:end,j,i))'], 3, 'g');
                mid_point = length(squeeze(K_w_m_yy_hist(nn,2:end,j,i)));
                K_e_sm_yy_hist(nn,1:size(K_e_m_yy_hist,2),j,i) = KK_hist((mid_point+1):end);
                K_w_sm_yy_hist(nn,2:size(K_w_m_yy_hist,2),j,i) = fliplr(KK_hist(1:mid_point));
            end
        end
    end
end

% Smooth RCP8.5 data
for nn = 1:num_years
    for j = 1:length(lat)
        for i = 1:size(K_e_m_yy_rcp,4)
            if size(K_w_m_yy_rcp,2) >= 3 && size(K_e_m_yy_rcp,2) >= 3
                KK_rcp = smooth([fliplr(squeeze(K_w_m_yy_rcp(nn,2:end,j,i))'), squeeze(K_e_m_yy_rcp(nn,1:end,j,i))'], 3, 'g');
                mid_point = length(squeeze(K_w_m_yy_rcp(nn,2:end,j,i)));
                K_e_sm_yy_rcp(nn,1:size(K_e_m_yy_rcp,2),j,i) = KK_rcp((mid_point+1):end);
                K_w_sm_yy_rcp(nn,2:size(K_w_m_yy_rcp,2),j,i) = fliplr(KK_rcp(1:mid_point));
            end
        end
    end
end

% Calculate zonal mean winds
uzm_hist = squeeze(mean(mean(uu_yy_hist(1:num_years,:,:,:),4),2));
uzm_rcp = squeeze(mean(mean(uu_yy_rcp(1:num_years,:,:,:),4),2));

save('uv_monthly_CESM2_hist_300hpa.mat','K_e_sm_yy_hist','K_w_sm_yy_hist','uzm_hist','lat');
save('uv_monthly_CESM2_rcp85_300hpa.mat','K_e_sm_yy_rcp','K_w_sm_yy_rcp','uzm_rcp','lat');

end

% Load and process historical data
load('uv_monthly_CESM2_hist_300hpa.mat','K_e_sm_yy_hist','K_w_sm_yy_hist','uzm_hist','lat');
uc_hist = squeeze(mean(uzm_hist,1))';
K_e_m_hist = squeeze(mean(K_e_sm_yy_hist,1));
K_w_m_hist = squeeze(mean(K_w_sm_yy_hist,1));

% Load and process RCP8.5 data
load('uv_monthly_CESM2_rcp85_300hpa.mat','K_e_sm_yy_rcp','K_w_sm_yy_rcp','uzm_rcp','lat');
uc_rcp85 = squeeze(mean(uzm_rcp,1))';
K_e_m_rcp = squeeze(mean(K_e_sm_yy_rcp,1));
K_w_m_rcp = squeeze(mean(K_w_sm_yy_rcp,1));

% Set up wavenumber and frequency arrays
max_wavenum = min(size(K_e_m_hist,3), size(K_e_m_rcp,3)) - 1;
wavenum = [0:max_wavenum]; 
freq = 0.5*[0:(size(K_w_m_hist,1)-1)]/(size(K_w_m_hist,1)-1); 
mm = 2:length(wavenum);
wavenum = wavenum(mm);
K_e_hist = K_e_m_hist(:,:,mm);
K_w_hist = K_w_m_hist(:,:,mm);
K_e_rcp = K_e_m_rcp(:,:,mm);
K_w_rcp = K_w_m_rcp(:,:,mm);

nk = length(wavenum);  
nj = length(lat);
kk = 1:nk; 
jj = 1:nj;
c_max = 50;  
c_min = -10;   
ca = 0.0:0.1:c_max;   
cc = 1:length(ca);

% Convert from frequency-wavenumber to phase speed-wavenumber space
for pos = 1:nj
    [K_e_new_hist, K_w_new_hist] = convert_fk_to_ck(freq, wavenum(kk), cos(lat(pos)/180.*pi), ca, K_e_hist(:,pos,kk), K_w_hist(:,pos,kk));
    K_e_lat_hist(:,pos) = sum(K_e_new_hist,2);
    K_w_lat_hist(:,pos) = sum(K_w_new_hist,2);
    
    [K_e_new_rcp, K_w_new_rcp] = convert_fk_to_ck(freq, wavenum(kk), cos(lat(pos)/180.*pi), ca, K_e_rcp(:,pos,kk), K_w_rcp(:,pos,kk));
    K_e_lat_rcp85(:,pos) = sum(K_e_new_rcp,2);
    K_w_lat_rcp85(:,pos) = sum(K_w_new_rcp,2);
end

save('co-spectra_cesm2_300hpa.mat','K_w_lat_hist','K_w_lat_rcp85','K_e_lat_hist','K_e_lat_rcp85','uc_hist','uc_rcp85');

% Plotting
figure('pos',[10 10 700 700]);
step = 0.01;
colorb = [-5*step:step:step*5];
[cc1,hh1] = contourf(-ca(cc),lat(jj),1.25*squeeze(K_w_lat_rcp85(cc,jj)-K_w_lat_hist(cc,jj))',colorb); hold on;
[cc2,hh2] = contourf(ca(cc),lat(jj),1.25*squeeze(K_e_lat_rcp85(cc,jj)-K_e_lat_hist(cc,jj))',colorb);
caxis([min(colorb),max(colorb)]);
gg = brewermap(10,'RdBu');
gg(5:6,:) = ones(2,3);
colormap(gg(end:-1:1,:)); 
set([hh1,hh2],'edgecolor','none');
step = 0.02;
[c1,h1] = contour(-ca(cc),lat(jj),squeeze(K_w_lat_hist(cc,jj))',step:step:step*12,'k-','linewidth',2); hold on;
[c2,h2] = contour(ca(cc), lat(jj),squeeze(K_e_lat_hist(cc,jj))',step:step:step*12,'k-','linewidth',2);
axis([c_min,c_max,-80,80]);
set([h1,h2],'linewidth',2.0)
xlabel('angular phase speed (m/s)','fontsize',20);
ylabel('latitude(deg)','fontsize',20);
set(gca,'xtick',[c_min:5:c_max],'ytick',[-80:20:80],'ylim',[-80,80],'fontsize',20);
line([0,0],[-80,80],'linestyle','-.','color','k');
plot(uc_rcp85./(cos(lat*pi/180)),lat,'r','linewidth',2);hold on;
plot(uc_hist./(cos(lat*pi/180)),lat,'k','linewidth',2);
title('Co-spectra Analysis at 300 hPa - Monthly Data','fontsize',16);