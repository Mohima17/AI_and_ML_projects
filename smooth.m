% S = smooth(F,nx,option);
% smooth the input function F with a window with a half-width nx
function S = smooth(F,nx,option);
if option == 'b'
  S = smooth_b(F,nx);
elseif option == 'g'
  S = smooth_g(F,nx);
end
return
 
%---------------------------------------------------------------
% box window
function B = smooth_b(F,nx)
len=length(F);
B=zeros(len,1);
for i=1:len
    B(i) = mean(F(max(1,i-nx):min(len,i+nx)));
end
 
%---------------------------------------------------------------
% Gaussian window 
function G = smooth_g(F,nx)
len=length(F);
G=zeros(len,1);
for i=1:len
for j=max(1,(i-3*nx)):min(len,(i+3*nx))
    G(i) = G(i)+ F(j)*exp(-0.5*((i-j)/nx)^2)/nx/sqrt(2*pi);
end
end