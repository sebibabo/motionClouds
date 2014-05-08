addpath('functions')
addpath(genpath('/home/bosse/dev/matlab/0.1-bosse/functionPool'))
clear all
path = 'mcForCmpWithMoco_coherence_highterSpatFreqOctBw_adap' ;


% take monitor parameters
% planar display
monitorRefresh = 60 ;
pixelPitch     = 0.2865/10 ;% pixelsize in cm

% experiments parametergenerateHighContrastStimuli_octSpBw
viewingDistance = 70 ;

frameRefresh   = monitorRefresh/2 ;
stimModRate    = 1 ;
numberOfReps   = 12 ;
framesPerMod   = frameRefresh/stimModRate/2 ; 

if ceil(framesPerMod) ~= framesPerMod
    disp('frames per modulation is not integer')
end


%%



vAnglePerSec    = 8.6;
vCoherent       = physicalSpeed2machineSpeed(vAnglePerSec,frameRefresh,viewingDistance,pixelPitch) ;
vIncoherent     = 0 ;


% only one fwhm
fwhm = [ 0.0338];



BvCoherent = fwhm/(2*sqrt(2*log(2)));
% BvIncoherent  = allBvCoherent(1);
BvIncoherent  = 10/(2*sqrt(2*log(2)));

%% set the the motion cloud parameters:
%  spatial properties

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

for coherence = [0.1 0.2 0.5 1] ;
  paramList{counter} = [ coherence ];
  counter = counter+1;
end

%%
zIncoherent     =  envelope_gabor_adap(fx, fy, ft, VxIncoherent, VyIncoherent,BvIncoherent, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
zCoherent =  envelope_gabor_adap(fx, fy, ft, VxCoherent, VyCoherent,BvCoherent, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;

zIncoherent = zIncoherent/sqrt(mean(zIncoherent(:).^2));
zCoherent   = zCoherent/sqrt(mean(zCoherent(:).^2));


for i = 1:length(paramList)
  
  
  coherence =  paramList{i}(1);
  
  filenamePrefix = sprintf('octBwRdkMc_%iHz_N_%i_coh_%.2f_Bsf_%.3f_sf_%.3f',frameRefresh ,Nx, coherence,B_sf,sf_0) ;
  % generate incoherent motion coud
  
  
  cloudIncoherent = rectif(random_cloud_adap(zIncoherent,1,0,0), 1,'Michelson',1) ;
  cloudIncoherent = cloudIncoherent(:,:,1:(framesPerMod-1)*numberOfReps);
  temp            = histeq(cloudIncoherent(:));
  cloudIncoherent = reshape( temp, size(cloudIncoherent ,1),size(cloudIncoherent ,2),size(cloudIncoherent ,3));
  
  sIncoherent = sum(zIncoherent,3) ; 
  
  %%
  zTotal          = sqrt( coherence*zCoherent.^2 + (1-coherence)*zIncoherent.^2) ;
    
  cloudTotal = rectif(random_cloud_adap(zTotal,1,0,0), 1,'Michelson',1) ;
  cloudTotal = cloudTotal(:,:,1:(framesPerMod+1)*numberOfReps);
  temp            = histeq(cloudTotal(:));
  cloudTotal = reshape( temp, size(cloudTotal,1),size(cloudTotal,2),size(cloudTotal,3));

  sTotal       = sum(zTotal,3) ; 
  
  %%
  % write out properties of clouds
  h = figure;
  subplot(4,2,1)
  imshow(cloudIncoherent(:,:,1)) 
  title(['1st frame of cloud, coherence = ',num2str(coherence)])
  
  subplot(4,2,2)
  imshow(cloudTotal(:,:,1)) 
  title(['1st frame of cloud, coherence = ',num2str(coherence)])
  
  subplot(4,2,3)
  hist(reshape(cloudIncoherent(:,:,1),1,[]),255)
  
  subplot(4,2,4)
  hist(reshape(cloudTotal(:,:,1),1,[]),255)
  
  subplot(4,2,5)
  plot(fx(1,:,1),sIncoherent(:,ceil(Nx/2)))
  hold on
  plot(fx(1,:,1),sTotal(:,ceil(Nx/2)),'r--')
  title('Projection to spatial frequencies of 1st frame')
%   legend({['[v, Bv] = ',num2str([VxIncoherent, BvIncoherent])],['[v, Bv] = ',num2str([VxCoherent, BvCoherent])]});

  subplot(4,2,6)
  plot(fx(1,:,1),sIncoherent(ceil(Nx/2),:))
  hold on
  plot(fx(1,:,1),sTotal(ceil(Nx/2),:),'r--')
  
  
  
  subplot(4,2,7)
  plot(fx(1,:,1),sIncoherent(:,ceil(Nx/2))./sTotal(:,ceil(Nx/2)))
  title('Ration between spatial frequencies ')

  subplot(4,2,8)
  plot(fx(1,:,1),sIncoherent(ceil(Nx/2),:)./sTotal(ceil(Nx/2),:))
    
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
%   cloudTotal(:,:,leftIdxs) = flipdim(cloudTotal(:,:,leftIdxs),3);   
%   
%   % up 
%   upIdxs = idxs(:,3:4:end);  % these are fine anyways
%   upIdxs = upIdxs (upIdxs ~=0);
%   cloudTotal(:,:,upIdxs ) = flipdim(permute(cloudTotal(:,:,upIdxs ),[2,1,3]),1);   
% 
% 
%   % down
%   downIdxs = idxs(:,4:4:end);  % these are fine anyways
%   downIdxs = downIdxs (downIdxs~=0);
%   cloudTotal(:,:,downIdxs ) = permute(cloudTotal(:,:,downIdxs ),[2,1,3]);   
  % change direction of motion in coherent part
  % right

  % left
  leftIdxs = idxs(:,2:4:end);  % these are fine anyways
  leftIdxs = leftIdxs(leftIdxs~=0);
  cloudTotal(:,:,leftIdxs) = flipdim(cloudTotal(:,:,leftIdxs),3);   
  
  % right
%   upIdxs = idxs(:,3:4:end);  % these are fine anyways
%   upIdxs = upIdxs (upIdxs ~=0);
%   cloudTotal(:,:,upIdxs ) = flipdim(permute(cloudTotal(:,:,upIdxs ),[2,1,3]),1);   


  % left
  downIdxs = idxs(:,4:4:end);  % these are fine anyways
  downIdxs = downIdxs (downIdxs~=0);
  cloudTotal(:,:,downIdxs ) = flipdim(cloudTotal(:,:,downIdxs ),3);     



% 
  
  
  
  mc_stim = zeros(Nx, Ny, framesPerMod*numberOfReps*2);
  idxs    = reshape(1:frameRefresh*numberOfReps,framesPerMod,[]);
  idxs_1  = idxs(:,1:2:end) ;
  idxs_2  = idxs(:,2:2:end) ;
  
  idxs_2=[idxs_1(end,:);idxs_2] ;
  idxs_1=idxs_1(1:end-1,:) ;
  
  mc_stim(:,:,idxs_1(:)) = cloudIncoherent ;
  mc_stim(:,:,idxs_2(:)) = cloudTotal ;
  mc_stim = drawFixationCross(mc_stim,30);
%   clear cloudIncoherent
  clear cloudTotal
  % put fixation cross
  
  
  images = single(reshape(mc_stim,size(mc_stim,1),size(mc_stim,2),1,size(mc_stim,3))) ;
  temp = 1:frameRefresh*numberOfReps ;
  temp = [temp; zeros(2,length(temp))];
  imageSequence = uint32(temp(:) );
        
  save(fullfile(path,[filenamePrefix, '.mat']),'images','imageSequence');
  write2aviByYuv(uint32(mc_stim*255), fullfile(path,[filenamePrefix, '.avi']), frameRefresh) ;
  
  
  write2aviByYuv(uint32(round(mc_stim)*255), fullfile(path,['binary_',filenamePrefix, '.avi']), frameRefresh) ;
end
%%

