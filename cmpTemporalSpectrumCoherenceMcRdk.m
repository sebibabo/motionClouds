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

VxCoherent = vx ;
VxIncoherent = 0 ;

VyCoherent   = 0 ;
VyIncoherent = 0 ;


%  spatial properties
sf_0 = 0.02;
B_sf = -1;%;0.001*6 ;

Nx = 512; 
Ny = Nx ;
N_frame = 45;% 2^size_T ;

alpha = 0 ;
theta = 0. ;
B_theta = inf ;
  loggabor = true ;
ft_0 = Nx/(N_frame) ;



%%
dataPath = '/export/data/motionClouds/MOCO_coh_MatFiles_30Hz_moreFrm' ;
allFilenames = {  'MOCO_coh000_F033_displ17_01.mat',
                  'MOCO_coh010_F033_displ17_02.mat',
                  'MOCO_coh020_F033_displ17_03.mat',
                  'MOCO_coh050_F033_displ17_04.mat',
                  'MOCO_coh100_F033_displ17_05.mat'} ;

% allFilenames = {  
%                   'MOCO_coh020_F033_displ17_03.mat',
%                   } ;                
%                 
                
% generate the two parts of the motion cloud
fwhm = [ 0.0338];
BvCoherent = fwhm/(2*sqrt(2*log(2)));
BvIncoherent  = 10/(2*sqrt(2*log(2)));                
[cloudsFx, cloudsFy, cloudsFt] = get_grids(Nx, Ny, N_frame)  ;

zIncoherent =  envelope_gabor_adap(cloudsFx, cloudsFy, cloudsFt, VxIncoherent, VyIncoherent,BvIncoherent, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
zCoherent   =  envelope_gabor_adap(cloudsFx, cloudsFy, cloudsFt, VxCoherent, VyCoherent,BvCoherent, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;

zIncoherent = zIncoherent/sqrt(mean(zIncoherent(:).^2));
zCoherent   = zCoherent/sqrt(mean(zCoherent(:).^2));

offset = 45 ;
thisLength = 15 ;
for fileIdx = 1:length(allFilenames)%1:1%
  h= figure('units','normalized','position',[0,0,1,1]); ;
  numberOfSpectra=4;
  %%
  % read the mocos
  % analyse bandwidth
  % plot spectra
  filename = allFilenames{fileIdx};
  load(fullfile(dataPath,filename))
  images =squeeze(images);

  part = 2;
  thisPartsImages = images(:,:,(part-1)*offset+1:(part-1)*offset+min(thisLength,offset));
  fx = ([1:size(thisPartsImages,1)]-size(thisPartsImages,1)/2)/size(thisPartsImages,1);
  ft = ([1:size(thisPartsImages,3)]-size(thisPartsImages,3)/2)/size(thisPartsImages,3);

  spectrumIdx = 1;
  spectrum        = abs(fftshift(fftn(thisPartsImages)));

  xySpectrum      = sum(spectrum,3);
  ytSpectrum      = squeeze(sum(spectrum,2)) ;
  xtSpectrum      = squeeze(sum(spectrum,1)) ;
  temporalFwhm    = getFwhm(xtSpectrum(length(fx)/2+1,:));

  if isnan(temporalFwhm)
    error('bv should not be nan!')
  end

  
  subplot(3,numberOfSpectra*2,(spectrumIdx-1)*2+1)
  imagesc(fx,fx,(xySpectrum))
  title('spatial spectrum as sum across ft, moco')
  xlabel('fx')
  ylabel('fy')

  subplot(3,numberOfSpectra*2,(spectrumIdx-1)*2+2)
  plot(fx,xySpectrum(length(fx)/2+1,:))
  grid on; xlim([-0.5,0.5]); ylim([0, max(xySpectrum(length(fx)/2+1,:))]);
  title('intersection along fx = 0');
  ylabel('fy');

  subplot(3,numberOfSpectra*2,numberOfSpectra*2+(spectrumIdx-1)*2+1)
  imagesc(fx,ft,ytSpectrum)
  title('spatiotemporal spectrum as sum across fy');
  xlabel('fx')
  ylabel('ft')

  subplot(3,numberOfSpectra*2,numberOfSpectra*2+(spectrumIdx-1)*2+2)
  plot(ft,ytSpectrum(length(fx)/2+1,:))
  grid on; xlim([-0.5,0.5]); ylim([0, max(ytSpectrum(length(fx)/2+1,:))]);
  title('intersection alog fx = 0');
  xlabel('ft')


  subplot(3,numberOfSpectra*2,numberOfSpectra*2*2+(spectrumIdx-1)*2+1)
  imagesc(fx,ft,xtSpectrum)
  title('spatiotemporal spectrum as sum across fx');
  xlabel('fy')
  ylabel('ft')

  subplot(3,numberOfSpectra*2,numberOfSpectra*2*2+(spectrumIdx-1)*2+2)
%   plot(ft,xtSpectrum(length(fx)/2+1,:))
  temp = fftshift(xtSpectrum);
  plot(ft,xtSpectrum(length(fx)/2+1,:))
% grid on; xlim([-0.5,0.5]); ylim([0, max(xtSpectrum(length(fx)/2+1,:))]);
%   plot(ft,temp(1,:))
  allMocoTemporalSpectra(fileIdx,:) = xtSpectrum(length(fx)/2+1,:);
   grid on; xlim([-0.5,0.5]); ylim([0, max(temp(1,:))]);
  title('intersection alog fy = 0');
  xlabel('ft')
  %%
  % generated cloud with this ratio of coherence
  % calculate coherence
  % options
  % a: take from filename, calc amp
%   coherence = str2num(filename(9:11))/100 ;
%   z         =  coherence*zCoherent + (1-coherence)*zIncoherent ;
  % b: take from filename, calc energy (THIS IS HOW IT WORKS!)
  coherence = str2num(filename(9:11))/100 ;
  z         = sqrt( coherence*zCoherent.^2 + (1-coherence)*zIncoherent.^2) ;
  % c: calc ratio of energies in moco
  % d: calc ratio of amplitudes in moco
  
     
  
  % plot spectra from freq domain
  
  

  fx = squeeze(cloudsFx(1,:,1));
  ft = squeeze(cloudsFt(1,1,:));
  % 
  %%
  spectrum = z;
  spectrumIdx = spectrumIdx +1;
  xySpectrum      = sum(spectrum,3);
  ytSpectrum      = squeeze(sum(spectrum,2)) ;
  xtSpectrum      = squeeze(sum(spectrum,1)) ;


  subplot(3,numberOfSpectra*2,(spectrumIdx-1)*2+1)
  imagesc(fx,fx,(xySpectrum))
  title('spatial spectrum as sum across ft, moco')
  xlabel('fx')
  ylabel('fy')

  subplot(3,numberOfSpectra*2,(spectrumIdx-1)*2+2)
  plot(fx,xySpectrum(length(fx)/2+1,:))
  grid on; xlim([-0.5,0.5]); ylim([0, max(xySpectrum(length(fx)/2+1,:))]);
  title('intersection along fx = 0');
  ylabel('fy');

  subplot(3,numberOfSpectra*2,numberOfSpectra*2+(spectrumIdx-1)*2+1)
  imagesc(fx,ft,ytSpectrum)
  title('spatiotemporal spectrum as sum across fy');
  xlabel('fx')
  ylabel('ft')

  subplot(3,numberOfSpectra*2,numberOfSpectra*2+(spectrumIdx-1)*2+2)
  plot(ft,ytSpectrum(length(fx)/2+1,:))
  grid on; xlim([-0.5,0.5]); ylim([0, max(ytSpectrum(length(fx)/2+1,:))]);
  title('intersection alog fx = 0');
  xlabel('ft')


  subplot(3,numberOfSpectra*2,numberOfSpectra*2*2+(spectrumIdx-1)*2+1)
  imagesc(fx,ft,xtSpectrum)
  title('spatiotemporal spectrum as sum across fx');
  xlabel('fy')
  ylabel('ft')

  subplot(3,numberOfSpectra*2,numberOfSpectra*2*2+(spectrumIdx-1)*2+2)
%   plot(ft,xtSpectrum(length(fx)/2+1,:))
%   grid on; xlim([-0.5,0.5]); ylim([0, max(xtSpectrum(length(fx)/2+1,:))]);
%   
  
  %   plot(ft,xtSpectrum(length(fx)/2+1,:))
%   temp = fftshift(xtSpectrum);
  plot(ft,xtSpectrum(length(fx)/2+1,:))
grid on; xlim([-0.5,0.5]); ylim([0, max(xtSpectrum(length(fx)/2+1,:))]);
%   plot(ft,temp(1,:))
  
  title('intersection alog fy = 0');
  xlabel('ft')

  % %%
  % transform to video domain
  % inv.trans to freq domin
  % plot spectra
  cloud           = random_cloud_adap(z,1,0,0);
  cloud           = cloud(:,:,1:thisLength);
  
  ft             = [1:thisLength]/thisLength-0.5;
  
  spectrum        = abs(fftshift(fftn(cloud)));

  spectrumIdx = spectrumIdx +1;
  xySpectrum      = sum(spectrum,3);
  ytSpectrum      = squeeze(sum(spectrum,2)) ;
  xtSpectrum      = squeeze(sum(spectrum,1)) ;
  
  subplot(3,numberOfSpectra*2,(spectrumIdx-1)*2+1)
  imagesc(fx,fx,(xySpectrum))
  title('spatial spectrum as sum across ft, moco')
  xlabel('fx')
  ylabel('fy')

  subplot(3,numberOfSpectra*2,(spectrumIdx-1)*2+2)
  plot(fx,xySpectrum(length(fx)/2+1,:))
  grid on; xlim([-0.5,0.5]); ylim([0, max(xySpectrum(length(fx)/2+1,:))]);
  title('intersection along fx = 0');
  ylabel('fy');

  subplot(3,numberOfSpectra*2,numberOfSpectra*2+(spectrumIdx-1)*2+1)
  imagesc(fx,ft,ytSpectrum)
  title('spatiotemporal spectrum as sum across fy');
  xlabel('fx')
  ylabel('ft')

  subplot(3,numberOfSpectra*2,numberOfSpectra*2+(spectrumIdx-1)*2+2)
  plot(ft,ytSpectrum(length(fx)/2+1,:))
  grid on; xlim([-0.5,0.5]); ylim([0, max(ytSpectrum(length(fx)/2+1,:))]);
  title('intersection alog fx = 0');
  xlabel('ft')


  subplot(3,numberOfSpectra*2,numberOfSpectra*2*2+(spectrumIdx-1)*2+1)
  imagesc(fx,ft,xtSpectrum)
  title('spatiotemporal spectrum as sum across fx');
  xlabel('fy')
  ylabel('ft')

  subplot(3,numberOfSpectra*2,numberOfSpectra*2*2+(spectrumIdx-1)*2+2)
  plot(ft,xtSpectrum(length(fx)/2+1,:))
  grid on; xlim([-0.5,0.5]); ylim([0, max(xtSpectrum(length(fx)/2+1,:))]);
  title('intersection alog fy = 0');
  xlabel('ft')

  % % % histeq contrast
  % % % inv.trans to freq domin
  % % % plot spectra
  rectifCloud           = rectif(cloud, 1,'Michelson',1) ;
  temp            = histeq(rectifCloud(:));
  histEqCloud          = reshape( temp, size(cloud,1),size(cloud,2),size(cloud,3));

  spectrum        = abs(fftshift(fftn(histEqCloud)));

  spectrumIdx = spectrumIdx +1;
  xySpectrum      = sum(spectrum,3);
  ytSpectrum      = squeeze(sum(spectrum,2)) ;
  xtSpectrum      = squeeze(sum(spectrum,1)) ;
  
  subplot(3,numberOfSpectra*2,(spectrumIdx-1)*2+1)
  imagesc(fx,fx,(xySpectrum))
  title('spatial spectrum as sum across ft, moco')
  xlabel('fx')
  ylabel('fy')

  subplot(3,numberOfSpectra*2,(spectrumIdx-1)*2+2)
  plot(fx,xySpectrum(length(fx)/2+1,:))
  grid on; xlim([-0.5,0.5]); ylim([0, max(xySpectrum(length(fx)/2+1,:))]);
  title('intersection along fx = 0');
  ylabel('fy');

  subplot(3,numberOfSpectra*2,numberOfSpectra*2+(spectrumIdx-1)*2+1)
  imagesc(fx,ft,ytSpectrum)
  title('spatiotemporal spectrum as sum across fy');
  xlabel('fx')
  ylabel('ft')

  subplot(3,numberOfSpectra*2,numberOfSpectra*2+(spectrumIdx-1)*2+2)
  plot(ft,ytSpectrum(length(fx)/2+1,:))
  grid on; xlim([-0.5,0.5]); ylim([0, max(ytSpectrum(length(fx)/2+1,:))]);
  title('intersection alog fx = 0');
  xlabel('ft')


  subplot(3,numberOfSpectra*2,numberOfSpectra*2*2+(spectrumIdx-1)*2+1)
  imagesc(fx,ft,xtSpectrum)
  title('spatiotemporal spectrum as sum across fx');
  xlabel('fy')
  ylabel('ft')

  subplot(3,numberOfSpectra*2,numberOfSpectra*2*2+(spectrumIdx-1)*2+2)
  temp = fftshift(xtSpectrum);
  plot(ft,xtSpectrum(length(fx)/2+1,:))
%   plot(ft,temp(1,:))
  grid on; xlim([-0.5,0.5]); ylim([0, max(xtSpectrum(length(fx)/2+1,:))]);
  title('intersection alog fy = 0');
  xlabel('ft')

  allMcTemporalSpectra(fileIdx,:) = xtSpectrum(length(fx)/2+1,:);
  plotFilename = sprintf('correctedCmpSpectraMcRdk_fwhm_%.3f_Bsf_%.3f_sf_%.3f',temporalFwhm, B_sf,sf_0) ;
  suptitle(filename);
  saveTightFigure(h, [plotFilename,'.pdf'])
%   close(h)
end
