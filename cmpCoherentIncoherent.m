% compare coherent and incoherent
addpath('functions')
addpath(genpath('/home/bosse/dev/matlab/0.1-bosse/functionPool'))
clear all
path = 'mcPilotStimuliCorrectedFwhm_highterSpatFreqOctBw_adap' ;

monitorRefresh = 60 ;
pixelPitch     = 0.2865/10 ;% pixelsize in cm

% experiments parameter
viewingDistance = 70 ;

frameRefresh   = monitorRefresh/2 ;
stimModRate    = 1 ;
numberOfReps   = 12 ;
framesPerMod   = frameRefresh/stimModRate/2 ; 



vAnglePerSec    = 8.6;
vCoherent       = physicalSpeed2machineSpeed(vAnglePerSec,frameRefresh,viewingDistance,pixelPitch) ;
vIncoherent     = 0 ;


allBvCoherent   = [0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5] ;
% fwhm from moco
fwhm = [10 0.6590 0.4425 0.3077 0.2712 0.2275 0.1839 0.0372];

BvCoherent = fwhm(end)/(2*sqrt(2*log(2)));
% BvIncoherent  = allBvCoherent(1);
BvIncoherent  = fwhm(1)/(2*sqrt(2*log(2)));
%% set the the motion cloud parameters:
%  spatial properties
sf_0 = 0.015;
B_sf = -1; % octave bandwidth


Nx = 512; 
Ny = Nx ;
N_frame = 180;% 2^size_T ;

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

% generate parameter list
paramList = {} ;
counter = 1;


%%
% generate incoherent motion coud
zIncoherent     =  envelope_gabor_adap(fx, fy, ft, VxIncoherent, VyIncoherent,BvIncoherent, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
cloudIncoherent = random_cloud_adap(zIncoherent,1,0,0) ;
cloudIncoherentConrastScaled = rectif(cloudIncoherent, 1,'Michelson',1) ;

temp            = histeq(cloudIncoherentConrastScaled (:));
cloudIncoherentHistEq = reshape( temp, size(cloudIncoherent ,1),size(cloudIncoherent ,2),size(cloudIncoherent ,3));
sIncoherent = sum(zIncoherent,3) ; 

% generate coherent motion cloud
zCoherent     =  envelope_gabor_adap(fx, fy, ft, VxCoherent, VyCoherent,BvCoherent, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
cloudCoherent = random_cloud_adap(zCoherent,1,0,0) ;
cloudCoherentConrastScaled = rectif(cloudCoherent, 1,'Michelson',1) ;

temp            = histeq(cloudCoherentConrastScaled (:));
cloudCoherentHistEq = reshape( temp, size(cloudCoherent ,1),size(cloudCoherent ,2),size(cloudCoherent ,3));
sCoherent = sum(zCoherent,3) ; 
%% the whole sequence
figure
subplot(3,2,1)
hist(cloudIncoherent(:),100)
subplot(3,2,2)
hist(cloudCoherent(:),100)

subplot(3,2,3)
hist(cloudIncoherentConrastScaled(:),100)
subplot(3,2,4)
hist(cloudCoherentConrastScaled(:),100)

subplot(3,2,5)
hist(cloudIncoherentHistEq(:),100)
subplot(3,2,6)
hist(cloudCoherentHistEq(:),100)
suplabel('along whole seq')
%% a random image
frameIdx = randi(180,1);
frameIncoherent = cloudIncoherent(:,:,frameIdx);
frameCoherent = cloudCoherent(:,:,frameIdx);

frameIncoherentConrastScaled = cloudIncoherentConrastScaled(:,:,frameIdx);
frameCoherentConrastScaled = cloudCoherentConrastScaled(:,:,frameIdx);

frameIncoherentHistEq = cloudIncoherentHistEq(:,:,frameIdx);
frameCoherentHistEq = cloudCoherentHistEq(:,:,frameIdx);


figure
subplot(3,2,1)
hist(frameIncoherent(:),100)
subplot(3,2,2)
hist(frameCoherent(:),100)

subplot(3,2,3)
hist(frameIncoherentConrastScaled(:),100)
subplot(3,2,4)
hist(frameCoherentConrastScaled(:),100)

subplot(3,2,5)
hist(frameIncoherentHistEq(:),100)
subplot(3,2,6)
hist(frameCoherentHistEq(:),100)
suplabel('rand imag')