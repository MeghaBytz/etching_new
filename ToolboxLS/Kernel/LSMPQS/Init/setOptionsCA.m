function options = setOptionsCA(varargin)

% setOptions
% options = setOptions('name1', value1, 'name2', value2, ...)
% options = setOptions(oldopts, 'name1', value1, 'name2', value2, ...)

% Creates an options structure (or alters an old one).
% Without any options, it displays the available options, values and
% defaults.


if (nargin == 0) && (nargout == 0)
  % Print out possible values of properties.
  fprintf('        accuracy: [ LSToolbox accuracy, medium by default ]\n');
  fprintf('               a: [ variable ]\n');
  fprintf('               b: [ variable ]\n');
  fprintf('              dc: [ curvature increment ]\n');
  
  fprintf('               f: [ factor in compressible model ]\n');
  fprintf('       Vmax_frac: [ target volume fraction in compressible model ]\n');
  
  fprintf('        GeomType: [ see setAvailGeom.m for available geometries]\n');
  fprintf('         ICplane: [ input condition interface ]\n');
  fprintf('          IC_pos: [ position of ''x'' or ''y'' IC plane ]\n');
  fprintf('              dx: [ grid spacing (in all dims) ]\n');
  
  fprintf('     segFileName: [ if GeomTYpe == segFile3DMA this has to be set to segmented file name]\n');
  fprintf('          magFac: [ magnification factor for segmented data]\n');
  fprintf('          doSeal: [ 1 if segmented data should be sealed on sides not relevant to flow]\n');  
  fprintf('          doFlip: [ 1 if (2D) geometry is to be flipped upside down (left-right) for ICplane = ''x''(''y'')]\n');  
  fprintf('            tMax: [ max running time allowed ]\n');
  fprintf('         epsStop: [ 1e-3 by default ]\n');
  fprintf('           tPlot: [ Time steps, relevant for plotting and/or oper. splitting ]\n');
  
  fprintf('        doReinit: [ 1 if reinitialization should be done after each step ]\n');  
  fprintf('     doOperSplit: [ 0 for no oper. splitting, 1 or 2 otherwise ]\n');
  
  fprintf('     OutFileName: [ output filename,  ''out'' by default\n');
  fprintf('          doSave: [ 1 if data, grid and mask should be saved ]\n');
  
  fprintf('       doDisplay: [ 1 if final results should be displayed ]\n');
  
  fprintf('        GeomType1: [ optional additional geometry]\n');
  fprintf('        GeomType2: [ optional additional geometry]\n');
  fprintf('        GeomType3: [ optional additional geometry]\n');
  fprintf('        PatchX   : [ 1 if GeomType1 to be patched in x direction, 2 if GeomType2 to be patched, 3 of both ]\n');
  fprintf('        PatchY   : [ 1 if GeomType1 to be patched in y direction, 2 if GeomType2 to be patched, 3 if both ]\n');
  fprintf('        addDuct  : [ 1 for adding a duct at the entrance, only x direction works for now]\n');
 
  fprintf('           doMask: [ 1 by default - masking out specified geometry]\n');
  fprintf('        stopTouch: [ 1 if the simulation is to be stopped if the opposite bdry is only touched.]\n');
  %fprintf('           doMin: [LSToolbox remnant, always set to 0]');
  
  fprintf('\n');
  return;
end

addpath(genpath('/ices/masha/ToolboxLS-1.0/Examples/Masha/'));
[AvailGeom, MaskGeom] = setAvailGeom;

%--------------------------------------------------------------------------
if((nargin > 0) && isstruct(varargin{1}))
    % First input argument is an old options structure
    options = varargin{1};
    startArg = 2;
else
    % Create the default options structure.
    options.accuracy = 'medium';
    options.a = 0.1;
    options.b = 0.05;
    options.dc = 0.1;
    options.ca = 0;
    options.velocity = 0;
    
    options.f = 1;
    options.Vmax_frac = 0.8;
    
    options.GeomType = AvailGeom(1).type;
    options.testGeom = 1;
    
    options.segFileName = '';
    options.magFac = 1;
    options.doSeal = 0;
    options.doFlip = 0;
    
    options.ICplane = 'x';
    options.IC_pos = 5;
    options.dx = 0.02;

    options.OutFileName = 'out';
    options.doSave = 1;

    options.epsStop = 1e-3;
    options.tMax = 15.0;
    options.tPlot = 0.1;

    options.doReinit = 1;
    options.doOperSplit = 0;
    options.doDisplay = 1;
    
    options.GeomType1 = AvailGeom(1).type;
    options.testGeom1 = 1;
    options.GeomType2 = AvailGeom(1).type;
    options.testGeom2 = 1;
    options.GeomType3 = AvailGeom(1).type;
    options.testGeom3 = 1;
    options.PatchX = 0;
    options.PatchY = 0;
    options.addDuct = 0;
    
    options.doMask = 1;
    options.doMin = 0;
    
    options.stopTouch = 1;

    startArg = 1;
end
  
%---------------------------------------------------------------------------
  % Loop through remaining name value pairs
for i = startArg : 2 : nargin
    name = varargin{i};
    value = varargin{i+1};

    % Remember that the case labels are lower case.
    switch(lower(name))
      case 'a'
        options.a = value;
      case 'b'
        options.b = value;
      case 'velocity'
        options.velocity = value;
      case 'ca'
        options.ca = value;
        
      case 'f'
      if(isa(value, 'double') && (prod(size(value)) == 1) && (value >= 0.0))
        options.f = value;
      else
        error('f must be a positive scalar double');
      end
      
      case 'vmax_frac'
      if(isa(value, 'double') && (prod(size(value)) == 1) && (value >= 0.0) && (value <= 1.0))
        options.Vmax_frac = value;
      else
        error('Vmax_frac must be a scalar double in [0,1]');
      end

     case 'da'
        options.da = value;
     
     case 'tmax'
      if(isa(value, 'double') && (prod(size(value)) == 1) && (value > 0.0))
        options.tMax = value;
      else
        error('tMax must be a positive scalar double value');
      end

      case 'epsstop'
      if(isa(value, 'double') && (prod(size(value)) == 1) && (value > 0.0))
        options.epsStop = value;
      else
        error('epsStop must be a positive scalar double value');
      end

      case 'tplot'
      if(isa(value, 'double') && (prod(size(value)) == 1) && (value > 0.0))
        options.tPlot = value;
      else
        error('tPlot must be a positive scalar double value');
      end

      case 'dx'
      if(isa(value, 'double') && (prod(size(value)) == 1) && (value > 0.0))
        options.dx = value;
      else
        error('dx must be a positive scalar double value');
      end

      case 'icplane'
      if( isa(value, 'char') && ( strcmp(value, 'x') || strcmp(value, 'y') || strcmp(value, 'z') ...
              || strcmp(value, 'c') || strcmp(value, 'b') || strcmp(value, 'r')) )
        options.ICplane = value;
      else
        error('ICplane must be one of the strings ''x'', ''y'', ''z'' or ''c''');
      end

      case 'ic_pos'
      %if( isa(value, 'int') )
        options.IC_pos = value;
      %else
      %  error('IC_pos must be a positive scalar int value');
      %end
      
      case 'outfilename'
      if isa(value, 'char')
          options.OutFileName = value;
      else
          error('OutFileName must be character array');
      end
      
      case 'segfilename'
      if isa(value, 'char')
          options.segFileName = value;
      else
          error('segFileName must be character array');
      end
      
      case 'magfac'
      if (value >=1 )
        options.magFac = value;
      else
        error('magFac must be positive integer');
      end
      
      case 'geomtype'  
      testGeom = 0; [m n] = size(AvailGeom);
      for i = 1:n
        if( strcmpi(value,AvailGeom(i).type) ) 
          options.testGeom = i;
          options.GeomType = AvailGeom(i).type;
        end
      end
      if (options.testGeom == 0)
        error('Uknown geometry type %s, reset to %s',value,AvailGeom(1).type);
        options.GeomType = AvailGeom(1).type;
        options.testGeom = 1;
      end
      
      case 'geomtype1'  
      testGeom1 = 0; [m n] = size(AvailGeom);
      for i = 1:n
        if( strcmpi(value,AvailGeom(i).type) ) 
          options.testGeom1 = i;
          options.GeomType1 = AvailGeom(i).type;
        end
      end
      
      if (options.testGeom1 == 0)
        error('Uknown geometry type %s, reset to %s',value,AvailGeom(1).type);
        options.GeomType1 = AvailGeom(1).type;
        options.testGeom1 = 1;
      end
      
      case 'geomtype2'  
      testGeom2 = 0; [m n] = size(AvailGeom);
      for i = 1:n
        if( strcmpi(value,AvailGeom(i).type) ) 
          options.testGeom2 = i;
          options.GeomType2 = AvailGeom(i).type;
        end
      end
      
      if (options.testGeom2 == 0)
        error('Uknown geometry type %s, reset to %s',value,AvailGeom(1).type);
        options.GeomType2 = AvailGeom(1).type;
        options.testGeom2 = 1;
      end
      
      case 'geomtype3'  
      testGeom3 = 0; [m n] = size(AvailGeom);
      for i = 1:n
        if( strcmpi(value,AvailGeom(i).type) ) 
          options.testGeom3 = i;
          options.GeomType3 = AvailGeom(i).type;
        end
      end
      
      if (options.testGeom3 == 0)
        error('Uknown geometry type %s, reset to %s',value,AvailGeom(1).type);
        options.GeomType3 = AvailGeom(1).type;
        options.testGeom3 = 1;
      end
      
      case 'patchx'
      if(isa(value, 'int') || (value == 0) || (value == 1) || (value == 2) || (value == 3) || (value == 4) )
        options.PatchX = value;
      else
        error('PatchX must be 0, 1 or 2');
      end
      
      case 'patchy'
      if(isa(value, 'int') || (value == 0) || (value == 1) || (value == 2) || (value == 3) || (value == 4))
        options.PatchY = value;
      else
        error('PatchY must be 0, 1 or 2');
      end
      
      case 'addduct'
      if(isa(value, 'int') || (value == 0) || (value == 1) )
        options.addDuct = value;
      else
        error('addDuct must be 0 or 1');
      end
      
      
     case 'doreinit'
      if(isa(value, 'int') || (value == 0) || (value == 1) )
        options.doReinit = value;
      else
        error('doReinit must be 0 or 1');
      end

      case 'doopersplit'
      if(isa(value, 'int') || (value == 0) || (value == 1) || (value == 2) )
        options.doOperSplit = value;
      else
        error('doOperSplit must be 0, 1 or 2');
      end

      case 'dosave'
      if(isa(value, 'int') || (value == 0) || (value == 1) )
        options.doSave = value;
      else
        error('doSave must be 0 or 1');
      end
      
      case 'doseal'
      if(isa(value, 'int') || (value == 0) || (value == 1) )
        options.doSeal = value;
      else
        error('doSeal must be 0 or 1');
      end
      
      case 'doflip'
      if(isa(value, 'int') || (value == 0) || (value == 1) )
        options.doFlip = value;
      else
        error('doFlip must be 0 or 1');
      end
      
      case 'dodisplay'
      if(isa(value, 'int') || (value == 0) || (value == 1) )
        options.doDisplay = value;
      else
        error('doDisplay must be 0 or 1');
      end
      
      case 'domask'
      if(isa(value, 'int') || (value == 0) || (value == 1) )
        options.doMask = value;
      else
        error('doMask must be 0 or 1');
      end
      
      case 'stoptouch'
      if(isa(value, 'int') || (value == 0) || (value == 1) )
        options.stopTouch = value;
      else
        error('stopTouch must be 0 or 1');
      end
      
      case 'accuracy'
       options.accuracy = 'high';    

     otherwise
      error('setOptions: Unknown option %s', name);

    end
end
