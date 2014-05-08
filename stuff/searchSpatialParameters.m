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

v              = 8 ;
framesPerMod   = frameRefresh/stimModRate/2 ; 


%% set the the motion cloud parameters:
% stimulus dependent 


% hardware dependent:
N_X = 512; 
N_Y = N_X ;
N_frame = N_X/4 ;


V_Y = 0. ;

alpha = 0 ;
theta = 0. ;
B_theta = inf ;
loggabor = true ;
ft_0 = 1;

if ceil(framesPerMod) ~= framesPerMod
    disp('frames per modulation is not integer')
end

[fx, fy, ft] = get_grids(N_X, N_Y, N_frame)  ;


% % generate parameter list
% paramList = {} ;
counter = 1;
Bv = 0.0001;
% for sf_0 = [ 0.005 0.01 0.05 0.1 ]
% for B_sf = [ 0.001, 0.01, 0.1  ]
%   paramList{counter} = [sf_0, B_sf];
%   counter = counter+1;
% end
% end

%sf_0 0.005 to coarse, 0.05 to fine
%B_sf 0.05 to noisy

  

for sf_0 = [ 0.015 0.02 ]
for B_sf = [ 0.001, 0.02  ]
  paramList{counter} = [sf_0, B_sf];
  counter = counter+1;
end
end




for paramsIdx = 1:length(paramList)
  
  sf_0 = paramList{paramsIdx}(1);
  B_sf = paramList{paramsIdx}(2);
  
  %   name = 'rdk_Vx_'+num2str(V_X)+'_Bv1_'+ num2str(Bv_1)+ '_Bv2_' + num2str(Bv)+ '_Bsf_' + num2str(Bsf);
  z{paramsIdx} =  envelope_gabor_adap(fx, fy, ft, v, V_Y,Bv, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
  cloud        = random_cloud(z{paramsIdx},1,0,0);
  rectifiedClouds{paramsIdx} = rectif(cloud, 1,'Michelson',1) ;
  
end
%%
for paramsIdx = 1:length(paramList)
  subplot(2,2,paramsIdx)
  imshow(rectifiedClouds{paramsIdx}(:,:,1)) 
  title(num2str(paramList{paramsIdx}))
%   implay(rectifiedClouds{paramsIdx})
%   title(num2str(paramList{paramsIdx}))
end
%%
%
% %%
% figure
% s1 = sum(z1,3) ;
% s2 = sum(z2,3) ;
% s3 = sum(z3,3) ;
% s4 = sum(z4,3) ;
% subplot(2,7,1)
% plot(s1(:,ceil(N_X/2)))
% hold on
% plot(s2(:,ceil(N_X/2)),'r')
% plot(s3(:,ceil(N_X/2)),'g')
% plot(s4(:,ceil(N_X/2)),'k')
% 
% subplot(2,7,2)
% imshow(rectifiedCloud1(:,:,1))
% 
% subplot(2,7,3)
% imshow(rectifiedCloud2(:,:,2))
% 
% subplot(2,7,4)
% imshow(rectifiedCloud3(:,:,3))
% 
% subplot(2,7,5)
% imshow(rectifiedCloud4(:,:,4))
% 
% subplot(2,7,6)
% d2 = s1./s2 ;
% d3 = s1./s3 ;
% d4 = s1./s4 ;
% plot(d2(:,ceil(N_X/2)),'r')
% hold on
% plot(d3(:,ceil(N_X/2)),'g')
% plot(d4(:,ceil(N_X/2)),'k')
% 
% subplot(2,7,7)
% plot(s1(ceil(N_X/2),:))
% hold on
% plot(s2(ceil(N_X/2),:),'r')
% plot(s3(ceil(N_X/2),:),'g')
% plot(s4(ceil(N_X/2),:),'k')
% 
% 
% rectZ1 = z1 ;
% rectZ2 = z2./repmat(sum(z2,3),[1 1 , size(z1,3)]).*repmat(sum(z1,3),[1 1 , size(z1,3)]);
% rectZ3 = z3./repmat(sum(z3,3),[1 1 , size(z1,3)]).*repmat(sum(z1,3),[1 1 , size(z1,3)]);
% rectZ4 = z4./repmat(sum(z4,3),[1 1 , size(z1,3)]).*repmat(sum(z1,3),[1 1 , size(z1,3)]);
% 
% s1 = sum(rectZ1,3) ;
% s2 = sum(rectZ2,3) ;
% s3 = sum(rectZ3,3) ;
% s4 = sum(rectZ4,3) ;
% subplot(2,7,8)
% plot(s1(:,ceil(N_X/2)))
% hold on
% plot(s2(:,ceil(N_X/2)),'r')
% plot(s3(:,ceil(N_X/2)),'g')
% plot(s4(:,ceil(N_X/2)),'k')
% 
% subplot(2,7,9)
% imshow(rectifiedRectCloud1(:,:,1))
% 
% subplot(2,7,10)
% imshow(rectifiedRectCloud2(:,:,2))
% 
% subplot(2,7,11)
% imshow(rectifiedRectCloud3(:,:,3))
% 
% subplot(2,7,12)
% imshow(rectifiedRectCloud4(:,:,4))
% 
% subplot(2,7,13)
% d2 = s1./s2 ;
% d3 = s1./s3 ;
% d4 = s1./s4 ;
% plot(d2(:,ceil(N_X/2)),'r')
% hold on
% plot(d3(:,ceil(N_X/2)),'g')
% plot(d4(:,ceil(N_X/2)),'k')
% 
% subplot(2,7,14)
% plot(s1(ceil(N_X/2),:))
% hold on
% plot(s2(ceil(N_X/2),:),'r')
% plot(s3(ceil(N_X/2),:),'g')
% plot(s4(ceil(N_X/2),:),'k')
% 
% 
