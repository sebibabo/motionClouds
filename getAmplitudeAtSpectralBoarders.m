addpath('functions')
addpath(genpath('/home/bosse/dev/matlab/0.1-bosse/functionPool'))
clear all

close all


% take monitor parameters
% planar display
monitorRefresh = 60 ;
pixelPitch     = 0.2865/10 ;% pixelsize in cm

% experiments parameter
viewingDistance = 70 ;

frameRefresh   = monitorRefresh/2 ;
stimModRate    = 1 ;
numberOfReps   = 12 ;
framesPerMod   = frameRefresh/stimModRate/2 ; 

if ceil(framesPerMod) ~= framesPerMod
    disp('frames per modulation is not integer')
end

vAnglePerSec    = 8.6;
vx              = physicalSpeed2machineSpeed(vAnglePerSec,frameRefresh,viewingDistance,pixelPitch) ;
vy              = 0;

%  spatial properties
sf_0 = 0.015*2;
B_sf = -1;%;0.001*6 ;

Nx = 512; 
Ny = Nx ;
N_frame = 180;% 2^size_T ;

alpha = 0 ;
theta = 0. ;
B_theta = inf ;
  loggabor = true ;
ft_0 = Nx/(N_frame) ;



%%
dataPath = '/export/data/motionClouds/MOCO_lifetime_MatFiles_30Hz_8.6deg_sec_MoreFrm' ;
allFilenames = {  'MOCO_random_F033_disp17_01.mat',
                  'MOCO_life002_F033_disp17_02.mat',
                  'MOCO_life003_F033_disp17_03.mat',
                  'MOCO_life004_F033_disp17_04.mat',
                  'MOCO_life005_F033_disp17_05.mat',
                  'MOCO_life006_F033_disp17_06.mat',
                  'MOCO_life007_F033_disp17_07.mat',
                  'MOCO_life008_F033_disp17_08.mat'};
                
allFilenames = {  'MOCO_random_F033_disp17_01.mat',
                  'MOCO_life008_F033_disp17_08.mat'};
                
allFilenames = {  
                  'MOCO_life008_F033_disp17_08.mat'};
% for sf_0 = 0.015*[1 1.5 1.6 1.7 1.8 1.93 1.95 2 ]
  for sf_0 = 0.015*[1.9 1.93 1.935 1.9353 1.94  ]
                
offset = 45 ;
for fileIdx = 1:length(allFilenames)%1:1%
  %%
  % read the mocos
  % analyse bandwidth
  filename = allFilenames{fileIdx};
  load(fullfile(dataPath,filename))
  images =squeeze(images);

  part = 2;
  thisPartsImages = images(:,:,(part-1)*offset+1:(part-1)*offset+offset);
  fx = ([1:size(thisPartsImages,1)]-size(thisPartsImages,1)/2)/size(thisPartsImages,1);
  ft = ([1:size(thisPartsImages,3)]-size(thisPartsImages,3)/2)/size(thisPartsImages,3);

  spectrumIdx = 1;
  spectrum        = abs(fftshift(fftn(thisPartsImages)));

  xySpectrum      = sum(spectrum,3);
  ytSpectrum      = squeeze(sum(spectrum,2)) ;
  xtSpectrum      = squeeze(sum(spectrum,1)) ;
  temporalFwhm    = getFwhm(xtSpectrum(length(fx)/2,:));

  if isnan(temporalFwhm)
    error('bv should not be nan!')
  end

  %%
  % generated cloud with this bandwidth
  % plot spectra from freq domain
  bv = temporalFwhm/(2*sqrt(2*log(2))); 

  [cloudsFx, cloudsFy, cloudsFt] = get_grids(Nx, Ny, N_frame)  ;
  z                              =  envelope_gabor_adap(cloudsFx, cloudsFy, cloudsFt, vx, vy,bv, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;

  fx = squeeze(cloudsFx(1,:,1));
  ft = squeeze(cloudsFt(1,1,:));
  % 
  %%
  spectrum = z;
  spectrumIdx = spectrumIdx +1;
  xySpectrum      = sum(spectrum,3);
  ytSpectrum      = squeeze(sum(spectrum,2)) ;
  xtSpectrum      = squeeze(sum(spectrum,1)) ;

  % get max. power ratio at the boarders
  maxPower = max(z(:).^2);
  powerRatio = z.^2/maxPower;
  maxPowerRatio(1) = max(max(max(powerRatio([1,end],:,:))));
  maxPowerRatio(2) = max(max(max(powerRatio(:,[1,end],:))));
  maxPowerRatio(3) = max(max(max(powerRatio(:,:,[1,end]))));
  
  maxAmplitude = max(z(:));
  amplitudeRatio = z/maxAmplitude;
  maxAmplitudeRatio(1) = max(max(max(amplitudeRatio([1,end],:,:))));
  maxAmplitudeRatio(2) = max(max(max(amplitudeRatio(:,[1,end],:))));
  maxAmplitudeRatio(3) = max(max(max(amplitudeRatio(:,:,[1,end]))));
  disp(['sf_0 = ', num2str(sf_0), ' max. power at boarder = ', num2str(max(maxPowerRatio)), ' max. amplitude at boarder = ', num2str(max(maxAmplitudeRatio))])
  
end
end