function [ yOut, schemeDataOut ] = maskAndKeepMin(t, yIn, schemeDataIn)
% maskAndKeepMin: Example postTimestep processing routine.
%
%  [ yOut, schemeDataOut ] = maskAndKeepMin(t, yIn, schemeDataIn)
%
%  This function demonstrates two entirely different processes that
%    can be accomplished through the postTimestep odeCFLn integrator option.
%
%  The first is to mask the evolving implicit surface function
%    (or otherwise modify its value after each timestep).
%
%  In this case masking is accomplished by
%
%          yOut = max(yIn, schemeDataIn.mask);
%
%
%   which ensures that phi cannot be negative anywhere that mask is positive.
%
%  The second is to keep track of some feature of that implicit surface
%    function or otherwise modify the schemeData structure after each
%    timestep.
%
%  In this case the feature recorded is the pointwise minimum over time of phi
%
%          schemeDataOut.min = min(yIn, schemeDataIn.min);
%
%
% Parameters:
%   t              Current time.
%   yIn            Input version of the level set function, in vector form.
%   schemeDataIn   Input version of a structure (see below).
%
%   yOut           Output version of the level set function, in vector form.
%   schemeDataOut  Output version of the structure (possibly modified).
%
% schemeData is a structure containing data specific to this type of 
%   term approximation.  For this function it contains the field(s)
%
%   .doMask      Boolean specifying whether masking should be performed.
%   .doMin       Boolean specifying whether min should be taken.
%   .mask	 Function against which to mask the level set function.
%   .min         Function which stores the minimum of the level set
%                  function over time (it is modified at each timestep).
%
% schemeData may contain other fields.

  checkStructureFields(schemeDataIn, 'doMask', 'doMin');

  % Mask the current level set function.
  if(schemeDataIn.doMask)
    checkStructureFields(schemeDataIn, 'mask');
    yOut = max(yIn, schemeDataIn.mask);
  else
    yOut = yIn;
  end

  % Record any new minimum values for each node.
  %   Use yOut to get the masked version of the data (if masking).
  schemeDataOut = schemeDataIn;
  if(schemeDataIn.doMin)
    checkStructureFields(schemeDataIn, 'min');
    schemeDataOut.min = min(yOut, schemeDataOut.min);
  end
