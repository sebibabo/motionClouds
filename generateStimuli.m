addpath('functions')
addpath(genpath('/home/bosse/dev/matlab/0.1-bosse/functionPool'))
clear all
path = 'mcPilotStimuli' ;


% take monitor parameters
% planar display
monitorRefresh = 60 ;
pixelPitch     = 0.2865/10 ;% pixelsize in cm

% experiments parameter
vAnglePerSec    = 4 ;
viewingDistance = 70 ;

frameRefresh   = monitorRefresh/3 ;
stimModRate    = 1 ;
numberOfReps   = 12 ;

v              = 3;%physicalSpeed2machineSpeed(vAnglePerSec,frameRefresh,viewingDistance,pixelPitch) ;
framesPerMod   = frameRefresh/stimModRate/2 ; 


%% set the the motion cloud parameters:
% stimulus dependent 
sf_0 = 0.05 ; % cycles per pix, should be changed to cycle per degree
B_sf  = 0.01;
Bv_1 = 1;
% hardware dependent:
N_X = 512; 
N_Y = N_X ;
N_frame = N_X;% 2^size_T ;
V_X = v ;
V_Y = 0. ;

alpha = 0 ;
theta = 0. ;
B_theta = inf ;
  loggabor = true ;
ft_0 = N_X/(N_frame) ;

if ceil(framesPerMod) ~= framesPerMod
    disp('frames per modulation is not integer')
end

[fx, fy, ft] = get_grids(N_X, N_Y, N_frame)  ;

if ~exist(path)
  mkdir(path)
end

% generate parameter list
paramList = {} ;
counter = 1;

for sf_0 = [0.01 0.1]
for Bv_2 = [ 0.001,  0.15, 0.3 ]
  
  paramList{counter} = [sf_0, Bv_2];
  counter = counter+1;
end
end

% for i = 1:length(paramList)
for i = 1:length(paramList)
  tic
  sf_0 = paramList{i}(1);
  Bv_2 = paramList{i}(2);
  
%   name = 'rdk_Vx_'+num2str(V_X)+'_Bv1_'+ num2str(Bv_1)+ '_Bv2_' + num2str(Bv_2)+ '_Bsf_' + num2str(Bsf);
  z1 =  envelope_gabor_adap(fx, fy, ft, V_X, V_Y,Bv_1, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
  cloud1 = rectif(random_cloud(z1,1,0,0), 1,'Michelson',1) ;
  cloud1 = cloud1(:,:,1:framesPerMod*numberOfReps);
  clear z1
  
  z2 =  envelope_gabor_adap(fx, fy, ft, V_X, V_Y,Bv_2, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
  cloud2 = rectif(random_cloud(z2,1,0,0), 1,'Michelson',1) ;
  cloud2 = cloud2(:,:,1:framesPerMod*numberOfReps);
%   
  clear z2
  %%
  idxs   = [1:framesPerMod*numberOfReps, zeros(1,mod(framesPerMod*numberOfReps,4))] ;
  idxs   = reshape(idxs,framesPerMod, []);
  % right
%   rightIdxs = idxs(:,1:4:end);  % these are fine anyways
%   rightIdxs = rightIdxs(rightIdxs~=0);
%   cloud2(:,:,rightIdxs) = cloud2(:,:,rightIdxs); 
  % left
  leftIdxs = idxs(:,2:4:end);  % these are fine anyways
  leftIdxs = leftIdxs(leftIdxs~=0);
  cloud2(:,:,leftIdxs) = flipdim(cloud2(:,:,leftIdxs),3);   
  
% up
  upIdxs = idxs(:,3:4:end);  % these are fine anyways
  upIdxs = upIdxs (upIdxs ~=0);
  cloud2(:,:,upIdxs ) = flipdim(permute(cloud2(:,:,upIdxs ),[2,1,3]),1);   


% down
  downIdxs = idxs(:,4:4:end);  % these are fine anyways
  downIdxs = downIdxs (downIdxs~=0);
  cloud2(:,:,downIdxs ) = permute(cloud2(:,:,downIdxs ),[2,1,3]);   
%%
  mc_stim = zeros(N_X, N_Y, framesPerMod*numberOfReps*2);
  idxs    = reshape(1:frameRefresh*numberOfReps,framesPerMod,[]);
  idxs_1  = idxs(:,1:2:end) ;
  idxs_2  = idxs(:,2:2:end) ;
  mc_stim(:,:,idxs_1(:)) = cloud1 ;
%   clear cloud1
  mc_stim(:,:,idxs_2(:)) = cloud2 ;
%   clear cloud2
  filenamePrefix = sprintf('rdkMc_N_%i_Bv1_%.1f_Bv2_%.3f_Bsf_%.3f_sf_%.3f',N_X, Bv_1,Bv_2,B_sf,sf_0) ;
  
  images = single(reshape(mc_stim,size(mc_stim,1),size(mc_stim,2),1,size(mc_stim,3))) ;
  temp = 1:frameRefresh*numberOfReps ;
  temp = [temp; zeros(2,length(temp))];
  imageSequence = uint32(temp(:) );
        
  save(fullfile(path,[filenamePrefix, '.mat']),'images','imageSequence');
  write2aviByYuv(uint32(mc_stim*255), fullfile(path,[filenamePrefix, '.avi']), frameRefresh) ;
  toc
end
