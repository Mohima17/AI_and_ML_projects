function dY = dy_dx(Y, x, dim)
% Calculate derivative dY/dx along dimension dim
% Y: input array
% x: coordinate vector
% dim: dimension along which to differentiate

if nargin < 3
    dim = 1;
end

% Get size of input array
sz = size(Y);
ndims_Y = ndims(Y);

% Permute array so that differentiation dimension comes first
if dim ~= 1
    perm_order = [dim, 1:dim-1, dim+1:ndims_Y];
    Y = permute(Y, perm_order);
    sz_perm = size(Y);
else
    sz_perm = sz;
    perm_order = 1:ndims_Y;
end

% Reshape for easier processing
n_along_dim = sz_perm(1);
n_other = prod(sz_perm(2:end));
Y_reshaped = reshape(Y, n_along_dim, n_other);

% Calculate derivatives
dY_reshaped = zeros(size(Y_reshaped));

% Use centered differences for interior points
for i = 2:n_along_dim-1
    dY_reshaped(i, :) = (Y_reshaped(i+1, :) - Y_reshaped(i-1, :)) / (x(i+1) - x(i-1));
end

% Use forward difference for first point
if n_along_dim > 1
    dY_reshaped(1, :) = (Y_reshaped(2, :) - Y_reshaped(1, :)) / (x(2) - x(1));
    % Use backward difference for last point
    dY_reshaped(end, :) = (Y_reshaped(end, :) - Y_reshaped(end-1, :)) / (x(end) - x(end-1));
end

% Reshape back
dY = reshape(dY_reshaped, sz_perm);

% Permute back to original dimension order
if dim ~= 1
    inv_perm_order(perm_order) = 1:ndims_Y;
    dY = permute(dY, inv_perm_order);
end

return