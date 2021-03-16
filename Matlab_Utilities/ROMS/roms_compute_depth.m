function [z_r, z_w, status] = roms_compute_depth(h,zeta,grid_info)
%=========================================================================%
% [z_r, z_w, status] = roms_compute_depth(h,zeta,grid_info)
% Computes depth (in m) at rho (r) and omega (w) points given bathymetric 
% depth 'h', sea surface height 'zeta' (can be 0), and the information in
% structure 'grid_info' which contains the following entries:
%  N           : # depth levels in rho grid
%  Vtransform  : Vertical transformation function index 
%  Vstretching : Vertical stretching function index
%  theta_s     : theta-surface parameter
%  theta_b     : theta-bottom parameter
%  hc          : aka TCline, aka critical depth
% 
% Note: h can be a single value, a vector, or 2d matrix.
% Note: zeta can either be a single value or an array the same size as h.
% Note: see ROMS wiki, etc., for help on what these variables mean in terms 
%       of ROMS' vertical coordinates
%=========================================================================%
% By Arin Nelson on 11/24/2020
% Last updated on 02/09/2021
%=========================================================================%

% Ensure grid_info contains the necessary variables
N           = grid_info.N;
Vtransform  = grid_info.Vtransform;
Vstretching = grid_info.Vstretching;
theta_s     = grid_info.theta_s;
theta_b     = grid_info.theta_b;
hc          = grid_info.hc;

% The sigma levels
k_w = 0:N;
k_r = (1:N)-0.5;
if(Vstretching<=4)
  s_w = (k_w-N) ./ N;
  s_r = (k_r-N) ./ N;
else
  s_w   = -(( k_w.^2 - 2*N.*k_w + k_w + N^2 - N)./(N^2 - N)) - 0.01.*( (k_w.^2 - N.*k_w) ./ (1-N) );
  s_r   = -(( k_r.^2 - 2*N.*k_r + k_r + N^2 - N)./(N^2 - N)) - 0.01.*( (k_r.^2 - N.*k_r) ./ (1-N) );
end

% Compute the vertical stretching function
switch Vstretching
  case 1
    Cs_r = (1-theta_b).*( sinh(theta_s.*s_r)./sinh(theta_s) ) + theta_b.*( (tanh(theta_s.*(s_r+0.5))./(2*tanh(0.5*theta_s))) - 0.5 ); 
    Cs_w = (1-theta_b).*( sinh(theta_s.*s_w)./sinh(theta_s) ) + theta_b.*( (tanh(theta_s.*(s_w+0.5))./(2*tanh(0.5*theta_s))) - 0.5 ); 
  case {2,3}
    gamma = 3;  % Set internally, see wiki page on Vertical_S-coordinate
    mu_r     = 0.5.*( 1 - tanh( gamma.*( s_r+0.5 ) ) );
    mu_w     = 0.5.*( 1 - tanh( gamma.*( s_w+0.5 ) ) );
    Cs_srf_r =  0 - log(cosh( gamma.*( abs(s_r).^theta_s ) ))./log(cosh(gamma));
    Cs_bot_r = -1 + log(cosh( gamma.*(    (s_r).^theta_b ) ))./log(cosh(gamma));
    Cs_srf_w =  0 - log(cosh( gamma.*( abs(s_w).^theta_s ) ))./log(cosh(gamma));
    Cs_bot_w = -1 + log(cosh( gamma.*(    (s_w).^theta_b ) ))./log(cosh(gamma));
    Cs_r     = mu_r.*Cs_bot_r + (1-mu_r).*Cs_srf_r;
    Cs_w     = mu_w.*Cs_bot_w + (1-mu_w).*Cs_srf_w;
  case {4,5}
    if(theta_s > 0)
      Cs_0_r = (1-cosh(theta_s.*s_r)) ./ (cosh(theta_s)-1);   
      Cs_0_w = (1-cosh(theta_s.*s_w)) ./ (cosh(theta_s)-1); 
    else
      Cs_0_r = -(s_r.^2);
      Cs_0_w = -(s_w.^2);
    end
    Cs_r = (exp(theta_b.*Cs_0_r) - 1) ./ (1-exp(-theta_b));
    Cs_w = (exp(theta_b.*Cs_0_w) - 1) ./ (1-exp(-theta_b));
  otherwise error(['Unknown vertical stretching function choice: ' num2str(Vstretching)]);
end

% Perform vertical transform for free-surface case
z_r = zeros(size(h,1),size(h,2),N);
z_w = zeros(size(h,1),size(h,2),N+1);
switch Vtransform
    case 1
      for i=1:N
        z_r(:,:,i) = hc*s_r(i) + (h-hc).*Cs_r(i);
      end
      for i=1:(N+1)
        z_w(:,:,i) = hc*s_w(i) + (h-hc).*Cs_w(i);
      end
    case 2
      for i=1:N
        z_r(:,:,i) = h.*( (hc*s_r(i) + h.*Cs_r(i)) ./ (hc + h) ); 
      end
      for i=1:(N+1)
        z_w(:,:,i) = h.*( (hc*s_w(i) + h.*Cs_w(i)) ./ (hc + h) ); 
      end
    otherwise
    error(['Possible valid Vtransform values are 1 and 2. Given: ' num2str(Vtransform)]);
end

% Squeezing out 1-length dimensions
z_r = squeeze(z_r);
z_w = squeeze(z_w);

% Return status=success
status = 1;

end