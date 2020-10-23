function decorr_time = r_compute_decorr_time(z,m,dt)
%=========================================================================%
% decorr_time = r_compute_decorr_time(z,m,dt)
% Computes the decorrelation time at each gridpoint in <z> in the same
% units as <dt>.  The mask <m> is used to skip land points.  If your ROMS
% didn't use masking, set m=ones(size(z,1),size(z,2)).
%=========================================================================%
% by Arin Nelson on 10/20/2020
%=========================================================================%

% Compute decorrelation time scales at each grid point
decorr_time = zeros(size(m));
for ix=1:size(m,1)
for iy=1:size(m,2)
if(m(ix,iy)==1)
  
  % Data at this grid point
  tmp = squeeze(z(ix,iy,:));
  if(~all(diff(tmp)==0) & ~all(isnan(tmp)))
 
    % For now, time series is 'detided' for testing purposes
    tmp = smooth(tmp,13,'lowess');

    % Autocorrelation
    acorr = autocorr(tmp,numel(tmp)-1);
  
    % Find first zero crossing
    ii = find(acorr(1:end-1)>0 & acorr(2:end)<=0,1,'first');
  
    % Save decorrelation time
    if(~isempty(ii))
      decorr_time(ix,iy) = dt*ii;
    else
      decorr_time(ix,iy) = 1e36;
    end
    
    % Clean-up
    clear acorr ii;
    
  end
  clear tmp;
    
end
end
end
clear ix iy;

% Debug: look at the map of decorr time's
%imagesc(decorr_time'); set(gca,'ydir','normal'); colorbar;

% DONE!
end