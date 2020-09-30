function in_vars = r_in_read(in_file)
%=========================================================================%
% in_vars = r_in_read(in_file)
% reads in the variable names and values from the roms.in file <in_file>
% and saves them to the structure in_vars, which is in the format of 
% in_vars(n).name  = < variable name >
% in_vars(n).value = < variable value(s) >
% 
% NOTE: For now, this assumes only 1 grid is being used.  This code will
% eventually be updated to accomodate multiple/nested grids.
%=========================================================================%
% By Arin Nelson on 09/23/2020
%=========================================================================%

% Open roms.in file
in_fid = fopen(in_file,'rt');
if(in_fid==-1)
  error(['Cannot find roms.in file: ' in_file]); 
end

% Read in lines of roms.in file and close it
in_line = {};
while(~feof(in_fid))
  in_line{end+1} = fgetl(in_fid);
end
fclose(in_fid);

% Remove comment lines (start with !) and blank lines
i = 1;
while(i<=numel(in_line))
  if(isempty(in_line{i}))
    in_line(i) = [];
  else
    if( strcmp(in_line{i}(1),'!') | all(isspace(in_line{i})) )
      in_line(i) = [];
    else
      i = i + 1;
    end
  end
end
clear i;

% Go through variables and get their values.
in_vars = struct;
for i=1:numel(in_line)
    
  % This line
  this_line = in_line{i};
  
  % Split via spaces
  this_split = strsplit(this_line,' ');
  
  % Remove empty entries
  n = 1;
  while(n<=numel(this_split))
    if(isempty(this_split{n}))
      this_split(n) = [];
    else
      n = n + 1;
    end
  end
  
  % First entry is the variable name
  in_vars(i).name = this_split{1};
  
  % If second entry is '=' or '==', remove it
  if(strcmp(this_split{2},'=') || strcmp(this_split{2},'=='))
    this_split{2} = [];
  end
  
  % The rest of the values are the variable values
  in_vars(i).value = this_split(3:end);

end %for i=1:numel(in_line)

end %function