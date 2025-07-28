% function [K_e_new, K_w_new] = convert_fk_to_ck(freq, wavenum, cos_lat, c, K_e, K_w)
% (freq, wavenum)::K_e, K_w --> (c, wavnum)::K_e_new, K_w_new
function [K_e_new, K_w_new] = convert_fk_to_ck(freq, wavenum, cos_lat, c, K_e, K_w)
RADIUS=6.371e6;          DAY=86400;           r_cos1=double(RADIUS*cos_lat);       
nk=length(wavenum);      nf=length(freq);     nc=length(c);
df=(2*pi)/(2*DAY)/(nf-1); 
[k_axis1,f_axis1]=meshgrid(wavenum,freq);
[k_axis2,c_axis2]=meshgrid(wavenum,c);
f_axis2 = (c_axis2*double(cos_lat)).*(k_axis2/r_cos1)/(2*pi/DAY);
% zonal wavenumber: k = 2*pi/L = n/r_cos
K_e_new=interp2(k_axis1,f_axis1,squeeze(K_e)/df,k_axis2,f_axis2).*(k_axis2/r_cos1)*double(cos_lat);
K_w_new=interp2(k_axis1,f_axis1,squeeze(K_w)/df,k_axis2,f_axis2).*(k_axis2/r_cos1)*double(cos_lat);
for k=1:nk
for cc=1:nc
  %dc=df*r_cos1/wavenum(k); % unresolved wave
  dc=df/wavenum(k);
  if(c(cc)<dc)
    K_e_new(cc,k)=nan;
    K_w_new(cc,k)=nan;
  else
    if(isnan(K_e_new(cc,k)))
      K_e_new(cc,k)=0;
    end
    if(isnan(K_w_new(cc,k)))
      K_w_new(cc,k)=0;
    end
  end  
end
end
return