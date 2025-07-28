% function f=dy_dx(y,x,n)
% dy/dx along n dimension
function f=dy_dx(y,x,n)
if(ndims(y) == 2)
   f=dy_dx_2d(y,x,n);
elseif(ndims(y) == 3)
   f=dy_dx_3d(y,x,n);
end
return
%----------------------------------------------------
function f=dy_dx_2d(y,x,n)
ni = size(y,n);
if (ni~=length(x))
  help dy_dx;
  disp('error: y and x should have equal length!');
end
 
if(n==1)
  dy_dx(1,:)=(y(1,:)-y(2,:))/(x(1)-x(2));
  dy_dx(ni,:)=(y(ni-1,:)-y(ni,:))/(x(ni-1)-x(ni));
  for i=2:(ni-1)
    dy_dx(i,:)=(y(i-1,:)-y(i+1,:))/(x(i-1)-x(i+1));
  end
elseif(n==2)
  dy_dx(:,1)=(y(:,1)-y(:,2))/(x(1)-x(2));
  dy_dx(:,ni)=(y(:,ni-1)-y(:,ni))/(x(ni-1)-x(ni));
  for i=2:(ni-1)
    dy_dx(:,i)=(y(:,i-1)-y(:,i+1))/(x(i-1)-x(i+1));
  end
end
f=dy_dx;
return
%----------------------------------------------------
function f=dy_dx_3d(y,x,n)
ni = size(y,n);
if (ni~=length(x))
  help dy_dx;
  disp('error: y and x should have equal length!');
end
 
if(n==1)
  dy_dx(1,:,:)=(y(1,:,:)-y(2,:,:))/(x(1)-x(2));
  dy_dx(ni,:,:)=(y(ni-1,:,:)-y(ni,:,:))/(x(ni-1)-x(ni));
  for i=2:(ni-1)
    dy_dx(i,:,:)=(y(i-1,:,:)-y(i+1,:,:))/(x(i-1)-x(i+1));
  end
elseif(n==2)
  dy_dx(:,1,:)=(y(:,1,:)-y(:,2,:))/(x(1)-x(2));
  dy_dx(:,ni,:)=(y(:,ni-1,:)-y(:,ni,:))/(x(ni-1)-x(ni));
  for i=2:(ni-1)
    dy_dx(:,i,:)=(y(:,i-1,:)-y(:,i+1,:))/(x(i-1)-x(i+1));
  end
elseif(n==3)
  dy_dx(:,:,1)=(y(:,:,1)-y(:,:,2))/(x(1)-x(2));
  dy_dx(:,:,ni)=(y(:,:,ni-1)-y(:,:,ni))/(x(ni-1)-x(ni));
  for i=2:(ni-1)
    dy_dx(:,:,i)=(y(:,:,i-1)-y(:,:,i+1))/(x(i-1)-x(i+1));
  end
end
f=dy_dx;
return