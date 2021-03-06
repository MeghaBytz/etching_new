function printOptions(options,fid)

% function printOptions(options,fid)
% Print out possible values in structure 'options' to file 'fid'
      
fprintf(fid,'\nstructure ''options''\n');
fprintf(fid,'   accuracy: %8s : [ LSToolbox accuracy, medium by default ]\n',options.accuracy);
fprintf(fid,'          a: %8g : [ constant; motion in normal direction ]\n',options.a);
fprintf(fid,'          b: %8g : [ constant; motion by mean curvature ]\n',options.b);
fprintf(fid,'         dc: %8g : [ curvature increment ]\n',options.dc);

fprintf(fid,'          f: %8g : [ factor in compressible model ]\n',options.f);
fprintf(fid,'  Vmax_frac: %8g : [ target volume fraction in compressible model ]\n',options.Vmax_frac);

fprintf(fid,'   GeomType: %8s : [ see setAvailGeom.m for available geometries ]\n',options.GeomType);
fprintf(fid,'segFileName: %8s : [if GeomTYpe == segFile3DMA has to be set to segmented file name]\n',options.segFileName);
fprintf(fid,'     magFac: %8d : [ magnification factor for segmented data]\n',options.magFac);
fprintf(fid,'     doSeal: %8d : [ 1 if segmented data should be sealed on sides not relevant to flow]\n',options.doSeal);
fprintf(fid,'     doFlip: %8d : [ 1 if (2D) geometry is to be flipped upside down (left-right) for ICplane = ''x''(''y'')]\n',options.doFlip);  

fprintf(fid,'  GeomType1: %8s : [ optional additional geometry]\n',options.GeomType1);
fprintf(fid,'  GeomType2: %8s : [ optional additional geometry]\n',options.GeomType2);
fprintf(fid,'  GeomType3: %8s : [ optional additional geometry]\n',options.GeomType3);
fprintf(fid,'     PatchX: %8d : [ 1 if Geomtype1 to be patched in x-dir, 2 if Geomtype2, 3 if both]\n',options.PatchX);
fprintf(fid,'     PatchY: %8d : [ 1 if Geomtype1 to be patched in y-dir, 2 if Geomtype2, 3 if both]\n',options.PatchY);
fprintf(fid,'    addDuct: %8d : [1 for adding a duct at the entrance, only x direction works for now]\n',options.addDuct);

if(options.ICplane == 'c' )
   fprintf(fid,'    ICplane: %8c : [ circular input condition interface ]\n',options.ICplane);
else
   fprintf(fid,'    ICplane: %6c=0 : [ input condition interface plane ]\n',options.ICplane);
end
fprintf(fid,'           IC_pos %d : [ position of ''x'' or ''y'' IC plane ]\n',options.IC_pos);
fprintf(fid,'         dx: %8g : [ grid spacing (in all dims) ]\n',options.dx);

fprintf(fid,'       tMax: %8g : [ max running time allowed ]\n',options.tMax);
fprintf(fid,'    epsStop: %8g : [ 1e-3 by default ]\n',options.epsStop);
fprintf(fid,'      tPlot: %8g : [ Time steps, relevant for plotting and/or oper. splitting ]\n',options.tPlot);

fprintf(fid,'   doReinit: %8d : [ 1 if reinitialization should be done after each step ]\n',options.doReinit);  
fprintf(fid,'doOperSplit: %8d : [ 0 for no oper. splitting, 1 or 2 otherwise ]\n',options.doOperSplit);

fprintf(fid,'OutFileName: %8s : [ output filename, dflt ''out'']\n',options.OutFileName);
fprintf(fid,'     doSave: %8d : [ 1 if data, grid and mask should be saved ]\n',options.doSave);

fprintf(fid,'  doDisplay: %8d : [ 1 if final results should be displayed ]\n',options.doDisplay);
fprintf(fid,'     doMask: %8d : [ masking out specified geometry; 1 by default]\n',options.doMask);
fprintf(fid,'     stopTouch: %8d : [stop simulation when the opposite side is touched]\n',options.stopTouch);
%fprintf(fid,'     doMin: %8d : [LSToolbox remnant, always set to 0]',options.doMin);

fprintf(fid,'\n');
  
