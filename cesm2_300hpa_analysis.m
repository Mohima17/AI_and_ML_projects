close all; clear; 
%param; 
% Modified for NetCDF files at 300 hPa level
 
RADIUS=6.371e6;
DAY=86400;
 
do_compute=1;  % Set to 1 to compute from NetCDF files
 
if do_compute>0 
%=========================== read data ================================
% Read historical data (1961-1980)
hist_file = 'b40.20th.track1.1deg.cam2.h0.1961-1980.1x1.nc';
% Read future data (2081-2100) 
rcp_file = 'b40.rcp8_5.1deg.cam2.h0.2081-2100.1x1.nc';

% Read coordinates first
lat = ncread(hist_file, 'lat');
lon = ncread(hist_file, 'lon');
lev = ncread(hist_file, 'lev');
time_hist = ncread(hist_file, 'time');

% Read coordinates from RCP file
lat_rcp = ncread(rcp_file, 'lat');
lon_rcp = ncread(rcp_file, 'lon');
lev_rcp = ncread(rcp_file, 'lev');
time_rcp = ncread(rcp_file, 'time');

% Level index for 300 hPa (index 16 as specified)
lev_idx = 16;

% Read wind data at 300 hPa level
% Historical period (1961-1980)
uu_hist = ncread(hist_file, 'U', [1, lev_idx, 1, 1], [Inf, 1, Inf, Inf]); % [lon, lat, time]
vv_hist = ncread(hist_file, 'V', [1, lev_idx, 1, 1], [Inf, 1, Inf, Inf]); % [lon, lat, time]

% Future period (2081-2100)
uu_rcp = ncread(rcp_file, 'U', [1, lev_idx, 1, 1], [Inf, 1, Inf, Inf]); % [lon, lat, time]
vv_rcp = ncread(rcp_file, 'V', [1, lev_idx, 1, 1], [Inf, 1, Inf, Inf]); % [lon, lat, time]

% Remove singleton dimension and permute to [time, lat, lon]
uu_hist = squeeze(uu_hist);  % [lon, lat, time]
vv_hist = squeeze(vv_hist);  % [lon, lat, time]
uu_rcp = squeeze(uu_rcp);    % [lon, lat, time]
vv_rcp = squeeze(vv_rcp);    % [lon, lat, time]

% Permute to [time, lat, lon] to match original code structure
uu_hist = permute(uu_hist, [3, 2, 1]);  % [time, lat, lon]
vv_hist = permute(vv_hist, [3, 2, 1]);  % [time, lat, lon]
uu_rcp = permute(uu_rcp, [3, 2, 1]);    % [time, lat, lon]
vv_rcp = permute(vv_rcp, [3, 2, 1]);    % [time, lat, lon]

rad = lat/180.*pi;

% Process historical data (20 years, assuming monthly data)
% Extract DJFM months for each year (Dec, Jan, Feb, Mar = 4 months per year)
nyears_hist = 20;  % 1961-1980
nyears_rcp = 20;   % 2081-2100

% Assuming monthly data (12 months per year)
for yy=1:nyears_hist
    % Extract DJFM months (months 12, 1, 2, 3) for each year
    if yy == 1
        % First year: only Jan, Feb, Mar (months 1, 2, 3)
        idx_start = 1;
        idx_end = 3;
    else
        % Other years: Dec from previous year + Jan, Feb, Mar from current year
        idx_start = (yy-2)*12 + 12;  % December from previous year
        idx_end = (yy-1)*12 + 3;     % March from current year
    end
    
    if idx_end <= size(uu_hist, 1)
        uu_yy_hist(yy,:,:,:) = uu_hist(idx_start:idx_end,:,:);
        vv_yy_hist(yy,:,:,:) = vv_hist(idx_start:idx_end,:,:);
    end
end

% Process RCP data similarly
for yy=1:nyears_rcp
    if yy == 1
        idx_start = 1;
        idx_end = 3;
    else
        idx_start = (yy-2)*12 + 12;
        idx_end = (yy-1)*12 + 3;
    end
    
    if idx_end <= size(uu_rcp, 1)
        uu_yy_rcp(yy,:,:,:) = uu_rcp(idx_start:idx_end,:,:);
        vv_yy_rcp(yy,:,:,:) = vv_rcp(idx_start:idx_end,:,:);
    end
end

%=========================== space-time analysis ================================
% Historical data processing
hann_win = hann(size(uu_yy_hist,2));  % Adjust window size to actual data length

for yy=1:nyears_hist   % beginning of loop for historical data
 
% ------------------ applying a hanning window ----------------------
for i=1:size(uu_yy_hist,4)
for j=1:size(uu_yy_hist,3)
for t=1:size(uu_yy_hist,2)
  u(i,t,j)=uu_yy_hist(yy,t,j,i)*hann_win(t);
  v(i,t,j)=vv_yy_hist(yy,t,j,i)*hann_win(t);
end
end
end
 
% ---------------------- space-time spectrum -------------------------
[East, West]=space_time_new(u,v);
 
% --------------------- momentum flux divergence ---------------------
for ff=1:size(East,1)
for jj=1:size(East,3)
K_e_mm_cos(ff,jj,:)=East(ff,:,jj)*(cos(rad(jj))^2);
K_w_mm_cos(ff,jj,:)=West(ff,:,jj)*(cos(rad(jj))^2);
end
end
dK_e_m=dy_dx(K_e_mm_cos,rad,2)/RADIUS*DAY;   % dy_dx: a subroutine for the meridional derivative
dK_w_m=dy_dx(K_w_mm_cos,rad,2)/RADIUS*DAY;
for ff=1:size(East,1)
for jj=1:size(East,3)
K_e_m_yy_hist(yy,ff,jj,:)=-dK_e_m(ff,jj,:)/(cos(rad(jj))^2);
K_w_m_yy_hist(yy,ff,jj,:)=-dK_w_m(ff,jj,:)/(cos(rad(jj))^2);
end
end
 
end           % end of loop for historical data

% RCP data processing
for yy=1:nyears_rcp   % beginning of loop for RCP data
 
% ------------------ applying a hanning window ----------------------
for i=1:size(uu_yy_rcp,4)
for j=1:size(uu_yy_rcp,3)
for t=1:size(uu_yy_rcp,2)
  u(i,t,j)=uu_yy_rcp(yy,t,j,i)*hann_win(t);
  v(i,t,j)=vv_yy_rcp(yy,t,j,i)*hann_win(t);
end
end
end
 
% ---------------------- space-time spectrum -------------------------
[East, West]=space_time_new(u,v);
 
% --------------------- momentum flux divergence ---------------------
for ff=1:size(East,1)
for jj=1:size(East,3)
K_e_mm_cos(ff,jj,:)=East(ff,:,jj)*(cos(rad(jj))^2);
K_w_mm_cos(ff,jj,:)=West(ff,:,jj)*(cos(rad(jj))^2);
end
end
dK_e_m=dy_dx(K_e_mm_cos,rad,2)/RADIUS*DAY;
dK_w_m=dy_dx(K_w_mm_cos,rad,2)/RADIUS*DAY;
for ff=1:size(East,1)
for jj=1:size(East,3)
K_e_m_yy_rcp(yy,ff,jj,:)=-dK_e_m(ff,jj,:)/(cos(rad(jj))^2);
K_w_m_yy_rcp(yy,ff,jj,:)=-dK_w_m(ff,jj,:)/(cos(rad(jj))^2);
end
end
 
end           % end of loop for RCP data
 
%=========================== smooth ================================
K_e_sm_yy_hist=zeros(size(K_e_m_yy_hist));
K_w_sm_yy_hist=zeros(size(K_w_m_yy_hist));
K_e_sm_yy_rcp=zeros(size(K_e_m_yy_rcp));
K_w_sm_yy_rcp=zeros(size(K_w_m_yy_rcp));

% removing stationary waves
K_e_m_yy_hist(:,1,:,:)=0;
K_w_m_yy_hist(:,1,:,:)=0;
K_e_m_yy_rcp(:,1,:,:)=0;
K_w_m_yy_rcp(:,1,:,:)=0;
 
% Smooth historical data
for nn=1:nyears_hist
for j=1:length(lat)
for i=1:size(K_e_m_yy_hist,4)
if size(K_e_m_yy_hist,2) > 2
    nfreq = size(K_e_m_yy_hist,2);
    KK=smooth([fliplr(squeeze(K_w_m_yy_hist(nn,2:nfreq,j,i))'),squeeze(K_e_m_yy_hist(nn,1:nfreq,j,i))'],3,'g');
    K_e_sm_yy_hist(nn,1:nfreq,j,i)=KK((nfreq):end);
    K_w_sm_yy_hist(nn,2:nfreq,j,i)=fliplr(KK(1:(nfreq-1)));
end
end
end
end

% Smooth RCP data
for nn=1:nyears_rcp
for j=1:length(lat)
for i=1:size(K_e_m_yy_rcp,4)
if size(K_e_m_yy_rcp,2) > 2
    nfreq = size(K_e_m_yy_rcp,2);
    KK=smooth([fliplr(squeeze(K_w_m_yy_rcp(nn,2:nfreq,j,i))'),squeeze(K_e_m_yy_rcp(nn,1:nfreq,j,i))'],3,'g');
    K_e_sm_yy_rcp(nn,1:nfreq,j,i)=KK((nfreq):end);
    K_w_sm_yy_rcp(nn,2:nfreq,j,i)=fliplr(KK(1:(nfreq-1)));
end
end
end
end

% Calculate mean zonal wind
uzm_hist=squeeze(mean(mean(uu_yy_hist,4),2));
uzm_rcp=squeeze(mean(mean(uu_yy_rcp,4),2));

% Save data
save('uv_DJFM_CESM2_hist_300hpa.mat','K_e_sm_yy_hist','K_w_sm_yy_hist','uzm_hist','lat');
save('uv_DJFM_CESM2_rcp85_300hpa.mat','K_e_sm_yy_rcp','K_w_sm_yy_rcp','uzm_rcp','lat');
 
end

% Load and process results
load('uv_DJFM_CESM2_hist_300hpa.mat','K_e_sm_yy_hist','K_w_sm_yy_hist','uzm_hist','lat');
uc_hist=squeeze(mean(uzm_hist,1))';
K_e_m_hist=squeeze(mean(K_e_sm_yy_hist,1));
K_w_m_hist=squeeze(mean(K_w_sm_yy_hist,1));
 
% Set up wavenumber and frequency arrays
nfreq = size(K_w_m_hist,1);
nwave = size(K_w_m_hist,3);
wavenum=[0:(nwave-1)]; 
freq=0.5*[0:(nfreq-1)]/(nfreq-1); 
mm=2:min(nwave,72);  % Adjust range based on actual data
wavenum = wavenum(mm);
K_e_hist = K_e_m_hist(:,:,mm);
K_w_hist = K_w_m_hist(:,:,mm);
nk=length(wavenum);  
nj=length(lat);
kk = 1:nk; 
jj=1:nj;
c_max=50;  
c_min = -10;   
ca=0.0:0.1:c_max;   
cc=1:length(ca);
 
% Convert frequency-wavenumber to phase speed-wavenumber for historical data
for pos=1:nj;
[K_e_new, K_w_new] = convert_fk_to_ck(freq, wavenum(kk), cos(lat(pos)/180.*pi), ca, K_e_hist(:,pos,kk), K_w_hist(:,pos,kk));
K_e_lat_hist(:,pos)=sum(K_e_new,2);
K_w_lat_hist(:,pos)=sum(K_w_new,2);
end
 
clear K_e_sm_yy_hist K_w_sm_yy_hist uzm_hist
 
% Process RCP data
load('uv_DJFM_CESM2_rcp85_300hpa.mat','K_e_sm_yy_rcp','K_w_sm_yy_rcp','uzm_rcp','lat');
uc_rcp85=squeeze(mean(uzm_rcp,1))';
K_e_m_rcp=squeeze(mean(K_e_sm_yy_rcp,1));
K_w_m_rcp=squeeze(mean(K_w_sm_yy_rcp,1));
K_e_rcp = K_e_m_rcp(:,:,mm);
K_w_rcp = K_w_m_rcp(:,:,mm);

% Convert for RCP data
for pos=1:nj;
[K_e_new, K_w_new] = convert_fk_to_ck(freq, wavenum(kk), cos(lat(pos)/180.*pi), ca, K_e_rcp(:,pos,kk), K_w_rcp(:,pos,kk));
K_e_lat_rcp85(:,pos)=sum(K_e_new,2);
K_w_lat_rcp85(:,pos)=sum(K_w_new,2);
end
 
save('co-spectra_cesm2_300hpa.mat','K_w_lat_hist','K_w_lat_rcp85','K_e_lat_hist','K_e_lat_rcp85','uc_hist','uc_rcp85');
 
% Plotting
figure('pos',[10 10 700 700]);
step = 0.01;
colorb=[-5*step:step:step*5];
[cc1,hh1]=contourf(-ca(cc),lat(jj),1.25*squeeze(K_w_lat_rcp85(cc,jj)-K_w_lat_hist(cc,jj))',colorb); hold on;
[cc2,hh2]=contourf(ca(cc),lat(jj),1.25*squeeze(K_e_lat_rcp85(cc,jj)-K_e_lat_hist(cc,jj))',colorb);
caxis([min(colorb),max(colorb)]);
gg=brewermap(10,'RdBu');
gg(5:6,:)=ones(2,3);
colormap(gg(end:-1:1,:)); 
set([hh1,hh2],'edgecolor','none');
step = 0.02;
[c1,h1]=contour(-ca(cc),lat(jj),squeeze(K_w_lat_hist(cc,jj))',step:step:step*12,'k-','linewidth',2); hold on;
[c2,h2]=contour(ca(cc), lat(jj),squeeze(K_e_lat_hist(cc,jj))',step:step:step*12,'k-','linewidth',2);
axis([c_min,c_max,-80,80]);
set([h1,h2],'linewidth',2.0)
xlabel('angular phase speed (m/s)','fontsize',20);
ylabel('latitude(deg)','fontsize',20);
title('Co-spectra Analysis at 300 hPa','fontsize',16);
set(gca,'xtick',[c_min:5:c_max],'ytick',[-80:20:80],'ylim',[-80,80],'fontsize',20);
line([0,0],[-80,80],'linestyle','-.','color','k');
plot(uc_rcp85./(cos(lat*pi/180)),lat,'r','linewidth',2);hold on;
plot(uc_hist./(cos(lat*pi/180)),lat,'k','linewidth',2);
legend('','','','','','RCP8.5 (2081-2100)','Historical (1961-1980)','Location','best');