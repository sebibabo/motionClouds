addpath('functions')
addpath(genpath('/home/bosse/dev/matlab/0.1-bosse/functionPool'))
clear all
path = 'mcPilotStimuliCorrectedFwhm_highterSpatFreqOctBw_adap' ;


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

vAnglePerSec    = 8.6;
vCoherent       = physicalSpeed2machineSpeed(vAnglePerSec,frameRefresh,viewingDistance,pixelPitch) ;
vIncoherent     = 0 ;


allBvCoherent   = [0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5] ;
% fwhm from moco
fwhm = [10 0.6590 0.4425 0.3077 0.2712 0.2275 0.1839 0.0372];
% fwhm = [10 0.6590 0.4425 0.3077 0.0372];
fwhm = [ 0.00372];
% fwhm = [0.6590 ];


allBvCoherent = fwhm/(2*sqrt(2*log(2)));
% BvIncoherent  = allBvCoherent(1);
BvIncoherent  = 10/(2*sqrt(2*log(2)));
%% set the the motion cloud parameters:
%  spatial properties
sf_0 = 0.015;
sf_0 = 0.02;
B_sf = -1; % octave bandwidth


Nx = 512; 
% moco: 18.5 deg
Ny = Nx ;
N_frame = 220;% 2^size_T ;

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


for BvCoherent = allBvCoherent
  paramList{counter} = [ BvCoherent ];
  counter = counter+1;
end

%%
% for i = 1:length(paramList)
for i = 1:length(paramList)
  tic
  
  BvCoherent  =  paramList{i}(1);
  
  filenamePrefix = sprintf('octBwRdkMc_%iHz_N_%i_Bv1_%.1f_Bv2_%.3f_Bsf_%.3f_sf_%.3f',frameRefresh ,Nx, BvIncoherent,BvCoherent,B_sf,sf_0) ;
  % generate incoherent motion coud
  zIncoherent     =  envelope_gabor_adap(fx, fy, ft, VxIncoherent, VyIncoherent,BvIncoherent, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
  cloudIncoherent = rectif(random_cloud_adap(zIncoherent,1,0,0), 1,'Michelson',1) ;
  cloudIncoherent = cloudIncoherent(:,:,1:(framesPerMod-1)*numberOfReps);
  temp            = histeq(cloudIncoherent(:));
  cloudIncoherent = reshape( temp, size(cloudIncoherent ,1),size(cloudIncoherent ,2),size(cloudIncoherent ,3));
  
  sIncoherent = sum(zIncoherent,3) ; 
%   clear zIncoherent
  % generate coherent motion cloud
  zCoherent =  envelope_gabor_adap(fx, fy, ft, VxCoherent, VyCoherent,BvCoherent, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
  %%
  
  sCoherent       = sum(zCoherent,3) ; 
  zCoherent=zCoherent./(repmat(sCoherent,[1,1,size(zCoherent,3)]));
  zCoherent=zCoherent.*(repmat(sIncoherent,[1,1,size(zCoherent,3)]));
  zCoherent(isnan(zCoherent)) = 0;
  %%
  cloudCoherent = rectif(random_cloud_adap(zCoherent,1,0,0), 1,'Michelson',1) ;
  cloudCoherent = cloudCoherent(:,:,1:(framesPerMod+1)*numberOfReps);

  temp            = histeq(cloudCoherent(:));
  cloudCoherent   = reshape( temp, size(cloudCoherent,1),size(cloudCoherent,2),size(cloudCoherent,3));
  
  sCoherent       = sum(zCoherent,3) ; 
  %%
  clear zCoherent
  % write out properties of clouds
  h = figure;
  subplot(4,2,1)
  imshow(cloudIncoherent(:,:,1)) 
  title(['1st frame of cloud, [v, Bv] = ',num2str([VxIncoherent, BvIncoherent])])
  
  subplot(4,2,2)
  imshow(cloudCoherent(:,:,1)) 
  title(['1st frame of cloud, [v, Bv] = ',num2str([VxCoherent, BvCoherent])])
  
  subplot(4,2,3)
  hist(reshape(cloudIncoherent(:,:,1),1,[]),255)
  
  subplot(4,2,4)
  hist(reshape(cloudCoherent(:,:,1),1,[]),255)
  
  subplot(4,2,5)
  plot(fx(1,:,1),sIncoherent(:,ceil(Nx/2)))
  hold on
  plot(fx(1,:,1),sCoherent(:,ceil(Nx/2)),'r--')
  title('Projection to spatial frequencies of 1st frame')
  legend({['[v, Bv] = ',num2str([VxIncoherent, BvIncoherent])],['[v, Bv] = ',num2str([VxCoherent, BvCoherent])]});

  subplot(4,2,6)
  plot(fx(1,:,1),sIncoherent(ceil(Nx/2),:))
  hold on
  plot(fx(1,:,1),sCoherent(ceil(Nx/2),:),'r--')
  
  
  
  subplot(4,2,7)
  plot(fx(1,:,1),sIncoherent(:,ceil(Nx/2))./sCoherent(:,ceil(Nx/2)))
  title('Ration between spatial frequencies ')

  subplot(4,2,8)
  plot(fx(1,:,1),sIncoherent(ceil(Nx/2),:)./sCoherent(ceil(Nx/2),:))
    
  saveTightFigure(h, fullfile(path,[filenamePrefix,'.pdf']));
  close all
  %%
  idxs   = [1:(framesPerMod+1)*numberOfReps, zeros(1,mod((framesPerMod+1)*numberOfReps,4))] ;
  idxs   = reshape(idxs,framesPerMod+1, []);
  
  % SB changing to right/left only
%   % change direction of motion in coherent part
%   % right
% 
%   % left
%   leftIdxs = idxs(:,2:4:end);  % these are fine anyways
%   leftIdxs = leftIdxs(leftIdxs~=0);
%   cloudCoherent(:,:,leftIdxs) = flipdim(cloudCoherent(:,:,leftIdxs),3);   
%   
%   % up 
%   upIdxs = idxs(:,3:4:end);  % these are fine anyways
%   upIdxs = upIdxs (upIdxs ~=0);
%   cloudCoherent(:,:,upIdxs ) = flipdim(permute(cloudCoherent(:,:,upIdxs ),[2,1,3]),1);   
% 
% 
%   % down
%   downIdxs = idxs(:,4:4:end);  % these are fine anyways
%   downIdxs = downIdxs (downIdxs~=0);
%   cloudCoherent(:,:,downIdxs ) = permute(cloudCoherent(:,:,downIdxs ),[2,1,3]);   
  % change direction of motion in coherent part
  % right

  % left
  leftIdxs = idxs(:,2:4:end);  % these are fine anyways
  leftIdxs = leftIdxs(leftIdxs~=0);
  cloudCoherent(:,:,leftIdxs) = flipdim(cloudCoherent(:,:,leftIdxs),3);   
  
  % right
%   upIdxs = idxs(:,3:4:end);  % these are fine anyways
%   upIdxs = upIdxs (upIdxs ~=0);
%   cloudCoherent(:,:,upIdxs ) = flipdim(permute(cloudCoherent(:,:,upIdxs ),[2,1,3]),1);   


  % left
  downIdxs = idxs(:,4:4:end);  % these are fine anyways
  downIdxs = downIdxs (downIdxs~=0);
  cloudCoherent(:,:,downIdxs ) = flipdim(cloudCoherent(:,:,downIdxs ),3);     



% 
  
  
  
  mc_stim = zeros(Nx, Ny, framesPerMod*numberOfReps*2);
  idxs    = reshape(1:frameRefresh*numberOfReps,framesPerMod,[]);
  idxs_1  = idxs(:,1:2:end) ;
  idxs_2  = idxs(:,2:2:end) ;
  
  idxs_2=[idxs_1(end,:);idxs_2] ;
  idxs_1=idxs_1(1:end-1,:) ;
  
  mc_stim(:,:,idxs_1(:)) = cloudIncoherent ;
  mc_stim(:,:,idxs_2(:)) = cloudCoherent ;
  mc_stim = drawFixationCross(mc_stim,30);
  clear cloudIncoherent
  clear cloudCoherent
  % put fixation cross
  
  
  images = single(reshape(mc_stim,size(mc_stim,1),size(mc_stim,2),1,size(mc_stim,3))) ;
  temp = 1:frameRefresh*numberOfReps ;
  temp = [temp; zeros(2,length(temp))];
  imageSequence = uint32(temp(:) );
        
  save(fullfile(path,[filenamePrefix, '.mat']),'images','imageSequence');
  write2aviByYuv(uint32(mc_stim*255), fullfile(path,[filenamePrefix, '.avi']), frameRefresh) ;
  toc
end
%%

