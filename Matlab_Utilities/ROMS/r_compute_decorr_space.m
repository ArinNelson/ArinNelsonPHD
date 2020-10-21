function decorr_space = r_compute_decorr_space(x,y,m,z)
%=========================================================================%
% results = r_compute_decorr_scales(x,y,m,z)
% Computes the decorrelation length scales for the variable <z> which has
% dimensions of (l,m,t) (meaning <x> and <y> are (l,m)) and land-sea mask
% <m>.  If your ROMS doesn't use masking, set m=ones(size(x)).
% The decorrelation length scales are computed by fitting temporal
% correlation vs. distance between grid points to a Gaussian in 4
% directions, NSEW, in that order, such that the output has dimensions
% (l,m,4).
%=========================================================================%
% by Arin Nelson on 10/20/2020
%=========================================================================%

% Possible user-changeable options in the future
min_to_calc = 3;     % Minimum number of grid points between data point and land to perform computation
max_to_calc = 100;   % Max number of grid points outwards to check for valid points to perform computation

% Compute decorrelation length scales at each grid point (N,E,S,W)
decorr_space = zeros(size(x,1),size(x,2),4);
for ix=1:size(x,1)
for iy=1:size(x,2)
clc; disp(['ix=' num2str(ix) ', iy=' num2str(iy)]);

  % Part 1: NORTHWARDS
  if(iy<size(x,2))
  if(sum(m(ix,iy:end))>=min_to_calc)
      
      % Gather data
      xN = x(ix,iy:end);
      yN = y(ix,iy:end);
      zN = squeeze(z(ix,iy:end,:));
      mN = m(ix,iy:end);
      
      % Only keep out to first land point
      iN = find(mN==0,1,'first');
      if(~isempty(iN))
        iN = iN-1;  xN = xN(1:iN);  yN = yN(1:iN);  zN = zN(1:iN,:);
      end
      clear iN mN;
      
      % Second check if enough data
      if(numel(xN)>=min_to_calc)
      
        % If too many, reduce
        if(numel(xN)>max_to_calc)
          xN = xN(1:max_to_calc); yN = yN(1:max_to_calc); zN = zN(1:max_to_calc,:);
        end
      
        % Compute distances
        dN = sqrt( (xN-xN(1)).^2 + (yN-yN(1)).^2 );
      
        % Compute correlations
        corrN = ones(numel(dN),1);
        for j=2:numel(dN)
          corrN(j) = nancorr(zN(1,:)',zN(j,:)');
        end
        clear j;
      
        % Compute fit to Gaussian correlation function, C(d)=exp(-d^2/2L^2)
        pN = polyfit(dN(:).^2,log(corrN(:)),1);
        lN = sqrt(-0.5/pN(1));
      
        % Save
        decorr_space(ix,iy,1) = lN;
        
        % Clean-up
        clear dN corrN pN lN;
      
      end
      clear xN yN zN;
      
  end
  end
  
  % Part 2: SOUTHWARDS
  if(iy>1)
  if(sum(m(ix,1:iy))>=min_to_calc)
      
      % Gather data
      xS = x(ix,1:iy);
      yS = y(ix,1:iy);
      zS = squeeze(z(ix,1:iy,:));
      mS = m(ix,1:iy);
      
      % Reverse order s.t. first point is data point at ix,iy
      xS = xS(end:-1:1);    yS = yS(end:-1:1);  mS = mS(end:-1:1);  zS = zS(end:-1:1,:);
      
      % Only keep out to first land point
      iS = find(mS==0,1,'first');
      if(~isempty(iS))
        iS = iS-1;  xS = xS(1:iS);  yS = yS(1:iS);  zS = zS(1:iS,:);
      end
      clear iS mS;
      
      % Second check if enough data
      if(numel(xS)>=min_to_calc)
      
        % If too many, reduce
        if(numel(xS)>max_to_calc)
          xS = xS(1:max_to_calc); yS = yS(1:max_to_calc); zS = zS(1:max_to_calc,:);
        end
      
        % Compute distances
        dS = sqrt( (xS-xS(1)).^2 + (yS-yS(1)).^2 );
      
        % Compute correlations
        corrS = ones(numel(dS),1);
        for j=2:numel(dS)
          corrS(j) = nancorr(zS(1,:)',zS(j,:)');
        end
        clear j;
      
        % Compute fit to Gaussian correlation function, C(d)=exp(-d^2/2L^2)
        pS = polyfit(dS(:).^2,log(corrS(:)),1);
        lS = sqrt(-0.5/pS(1));
      
        % Save
        decorr_space(ix,iy,2) = lS;
      
        % Clean-up
        clear dS corrS pS lS;
      
      end
      clear xS yS zS;
      
  end
  end
  
  % Part 3: EASTWARDS
  if(ix<size(x,1))
  if(sum(m(ix:end,iy))>=min_to_calc)
    
      % Gather data
      xE = x(ix:end,iy);
      yE = y(ix:end,iy);
      zE = squeeze(z(ix:end,iy,:));
      mE = m(ix:end,iy);
      
      % Only keep out to first land point
      iE = find(mE==0,1,'first');
      if(~isempty(iE))
        iE = iE-1;  xE = xE(1:iE);  yE = yE(1:iE);  zE = zE(1:iE,:);
      end
      clear iE mE;
      
      % Second check if enough data
      if(numel(xE)>=min_to_calc)     
      
        % If too many, reduce
        if(numel(xE)>max_to_calc)
          xE = xE(1:max_to_calc); yE = yE(1:max_to_calc); zE = zE(1:max_to_calc,:);
        end
      
        % Compute distances
        dE = sqrt( (xE-xE(1)).^2 + (yE-yE(1)).^2 );
      
        % Compute correlations
        corrE = ones(numel(dE),1);
        for j=2:numel(dE)
          corrE(j) = nancorr(zE(1,:)',zE(j,:)');
        end
        clear j;
      
        % Compute fit to Gaussian correlation function, C(d)=exp(-d^2/2L^2)
        pE = polyfit(dE(:).^2,log(corrE(:)),1);
        lE = sqrt(-0.5/pE(1));
      
        % Save
        decorr_space(ix,iy,3) = lE;
      
        % Clean-up
        clear dE corrE pE lE;
      
      end
      clear xE yE zE;
        
  end
  end
  
  % Part 4: WESTWARDS
  if(iy<size(x,2))
  if(sum(m(1:ix,iy))>=min_to_calc)
      
      % Gather data
      xW = x(1:ix,iy);      yW = y(1:ix,iy);    mW = m(1:ix,iy);    zW = squeeze(z(1:ix,iy,:));
      
      % Reverse order s.t. first point is data point at ix,iy
      xW = xW(end:-1:1);    yW = yW(end:-1:1);  mW = mW(end:-1:1);  zW = zW(end:-1:1,:);
      
      % Only keep out to first land point
      iW = find(mW==0,1,'first');
      if(~isempty(iW))
        iW = iW-1;  xW = xW(1:iW);  yW = yW(1:iW);  zW = zW(1:iW,:);
      end
      clear iW mW;
      
      % Second check if enough data
      if(numel(xW)>=min_to_calc)
      
        % If too many, reduce
        if(numel(xW)>max_to_calc)
          xW = xW(1:max_to_calc); yW = yW(1:max_to_calc); zW = zW(1:max_to_calc,:);
        end
      
        % Compute distances
        dW = sqrt( (xW-xW(1)).^2 + (yW-yW(1)).^2 );
      
        % Compute correlations
        corrW = ones(numel(dW),1);
        for j=2:numel(dW)
          corrW(j) = nancorr(zW(1,:)',zW(j,:)');
        end
        clear j;
      
        % Compute fit to Gaussian correlation function, C(d)=exp(-d^2/2L^2)
        pW = polyfit(dW(:).^2,log(corrW(:)),1);
        lW = sqrt(-0.5/pW(1));
      
        % Save
        decorr_space(ix,iy,4) = lW;
      
        % Clean-up
        clear dW corrW pW lW;
        
      end
      clear xW yW zW
      
  end
  end
    
end
end
clear ix iy;

% Debug check: plots?
%for i=1:4; subplot(2,2,i); imagesc(real(decorr_space(:,:,i)')); set(gca,'ydir','normal'); colorbar; caxis([0 1e-8]); end
%imagesc(real(sum(decorr_space,3)./4)'); set(gca,'ydir','normal'); colorbar; caxis([0 1e-8]);

% DONE!
end