function Results = r_nc_sanitycheck(fileIn)
% Results = r_nc_sanitycheck(fileName)
% Reads through the variables in the NetCDF file <fileName> and looks for 
% missing (NaN) values.

% Check that netcdf file exists
try info = ncinfo(fileIn);
catch err; error(['Could not find NetCDF file fileIn: ' fileIn]);
end

% Init the Results container
Results = cell(numel(info.Variables),1);

% If it does exist, begin looping through the variables
for i=1:numel(info.Variables)
    
  % If variable doesn't have any dimensions, it's a single value
  if(isempty(info.Variables(i).Dimensions))
    
    % Read in the value.  If it's NaN, replace it with fillValue
    varValue = ncread(fileIn,info.Variables(i).Name);
    if(isnan(varValue))
      Results{i} = ['Single-valued variable ' info.Variables(i).Name ' missing!.'];
    else
      Results{i} = 'ok';
    end
    clear varValue;
      
  else   
      
    % Get info on variable dimensions
    varDims = {info.Variables(i).Dimensions.Name};
    iDim    = zeros(numel(varDims),1);
    for j=1:numel(varDims)
      iDim = find(strcmp(varDims(j),{info.Variables(i).Dimensions.Name})==1);
    end
    
    % If no dimensions are 'time' or 'unlimited', read in the full variable
    if(~any([info.Dimensions(iDim).Unlimited]==1))
        
      % Read in the variable  If any values are NaN, report on it.
      varValue = ncread(fileIn,info.Variables(i).Name);
      if(isnan(varValue))
        Results{i} = ['Variable ' info.Variables(i).Name ' has missing values!'];
      else
        Results{i} = 'ok';
      end
      clear varValue;
        
    else
        
      % Loop through the unlimited dimension (always the last one)
      for j=1:info.Variables(i).Dimensions(end).Length
        
        % Dimension array
        dimArray = info.Variables(i).Size;
        dimArray(end) = j;
        
        % Load this variable
        test = ncread(fileIn,info.Variables(i).Name,dimArray);
        
        % Check if any are NaN.  If so, report on it.
        if(any(isnan(test)))
          Results{i} = [Results{i} num2str(j) ', '];
        end
        clear dimArray test;
          
      end
      clear j;
      
      % If no missing values found, say so.
      if(isempty(Results{i}))
        Results{i} = 'ok';
      end
        
    end
        
  end 
    
end
clear i;

end