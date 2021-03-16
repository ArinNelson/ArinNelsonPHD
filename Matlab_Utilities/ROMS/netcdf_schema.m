function schema = netcdf_schema

  % The info structure
  schema = struct;

  % Default information
  schema.Filename = '';     % Blank, will be set by some other code
  schema.Name '/';          % Default
  schema.Groups = [];       % Default
  schema.Format = '64bit';  % Default

  % Internal structures
  schema.Dimensions = struct;
  schema.Variables  = struct;
  schema.Attributes = struct;

end