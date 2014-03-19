function [sf0 Bsf, v,Bv] = physicalUnits2discreteUnits(sf0_cpd, Bsf_cpd, v_dps, B_dps,pixelPitch,viewingDistance,frameRate)
%   % laptop monitor
%   pixelPitch      = 0.22/10; % in cm
%   viewingDistance = 50;
%   sf0_cpd         = 4 ;
%   Bsf_cpd         = 1 ;
%   v_dps           = 4 ;
%   B_dps           = 1 ;
%   frameRate       = 20 ;
  % convert to machine units
  cmPerDegree    = 2*viewingDistance*tand(1/2);
  pxPerDegree    = cmPerDegree/pixelPitch;

  sf0            = sf0_cpd/pxPerDegree;
  Bsf            = Bsf_cpd/pxPerDegree;

  v              = v_dps/frameRate*pxPerDegree;
  Bv             = B_dps/frameRate*pxPerDegree;
