addpath('functions')
addpath(genpath('/home/bosse/dev/matlab/0.1-bosse/functionPool'))
clear all
path = 'mcPilotStimuli' ;


% take monitor parameters
% planar display
monitorRefresh = 60 ;
pixelPitch     = 0.2865/10 ;% pixelsize in cm

% experiments parameter
viewingDistance = 70 ;

frameRefresh   = monitorRefresh/3 ;
stimModRate    = 1 ;
numberOfReps   = 12 ;
framesPerMod   = frameRefresh/stimModRate/2 ; 

if ceil(framesPerMod) ~= framesPerMod
    disp('frames per modulation is not integer')
end

vAnglePerSec    = 4 ;
vCoherent       = 8;%physicalSpeed2machineSpeed(vAnglePerSec,frameRefresh,viewingDistance,pixelPitch) ;
vIncoherent     = 0 ;

BvIncoherent    = 0.5 ;
allBvCoherent   = [0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5] ;

%% set the the motion cloud parameters:
%  spatial properties
sf_0 = 0.015;
B_sf = 0.001 ;

Nx = 1024; 
Ny = Nx ;
N_frame = 128;% 2^size_T ;

VxCoherent   = vCoherent ;
VxIncoherent = vIncoherent;
VyCoherent   = 0. ;
VyIncoherent   = 0. ;

alpha = 0 ;
theta = 0. ;
B_theta = inf ;
  loggabor = true ;
ft_0 = Nx/(N_frame) ;

[fx, fy, ft] = get_grids(Nx, Ny, N_frame)  ;

if ~exist(path)
  mkdir(path)
end

  BvCoherent  =  0.0001;
  
  
  % generate coherent motion cloud
  zCoherent =  envelope_gabor_adap(fx, fy, ft, VxCoherent, VyCoherent,BvCoherent, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
  cloudCoherent = rectif(random_cloud(zCoherent,1,0,0), 1,'Michelson',1) ;
