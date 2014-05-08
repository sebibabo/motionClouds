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

frameRefresh   = monitorRefresh/2 ;
stimModRate    = 1 ;
numberOfReps   = 12 ;

v1              = physicalSpeed2machineSpeed(0,frameRefresh,viewingDistance,pixelPitch) ;
v2              = physicalSpeed2machineSpeed(1,frameRefresh,viewingDistance,pixelPitch) ;
v3              = physicalSpeed2machineSpeed(4,frameRefresh,viewingDistance,pixelPitch) ;
v4              = physicalSpeed2machineSpeed(10,frameRefresh,viewingDistance,pixelPitch) ;
framesPerMod     = frameRefresh/stimModRate/2 ; 


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
% counter = 1;
% 
% for sf_0 = [ 0.01 0.011 0.012 0.013 0.014]
% for Bv = [ 0.001]%,  0.15, 0.3 ]
%   
%   paramList{counter} = [sf_0, Bv];
%   counter = counter+1;
% end
% end


tic
sf_0  = 0.01 ;
B_sf  = 0.001;
Bv    = 0.11 ;


%   name = 'rdk_Vx_'+num2str(V_X)+'_Bv1_'+ num2str(Bv_1)+ '_Bv2_' + num2str(Bv)+ '_Bsf_' + num2str(Bsf);
z1 =  envelope_gabor_adap(fx, fy, ft, v1, V_Y,Bv, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
cloud1 = random_cloud(z1,1,0,0);
rectifiedCloud1 = rectif(cloud1, 1,'Michelson',1) ;

z2 =   envelope_gabor_adap(fx, fy, ft, v2, V_Y,Bv, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
cloud2 = random_cloud(z2,1,0,0) ;
rectifiedCloud2 = rectif(cloud2, 1,'Michelson',1) ;

z3 =  envelope_gabor_adap(fx, fy, ft, v3, V_Y,Bv, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
cloud3 = random_cloud(z3,1,0,0) ;
rectifiedCloud3 = rectif(cloud3, 1,'Michelson',1) ;

z4 =  envelope_gabor_adap(fx, fy, ft, v4, V_Y,Bv, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
cloud4 = random_cloud(z4,1,0,0) ;
rectifiedCloud4 = rectif(cloud4, 1,'Michelson',1) ;


rectZ1 = z1 ;
rectZ2 = z2./repmat(sum(z2,3),[1 1 , size(z1,3)]).*repmat(sum(z1,3),[1 1 , size(z1,3)]);
rectZ3 = z3./repmat(sum(z3,3),[1 1 , size(z1,3)]).*repmat(sum(z1,3),[1 1 , size(z1,3)]);
rectZ4 = z4./repmat(sum(z4,3),[1 1 , size(z1,3)]).*repmat(sum(z1,3),[1 1 , size(z1,3)]);

rectZ1(isnan(rectZ1)) = 0;
rectZ2(isnan(rectZ2)) = 0;
rectZ3(isnan(rectZ3)) = 0;
rectZ4(isnan(rectZ4)) = 0;

rectCloud1 = random_cloud(rectZ1,1,0,0);
rectifiedRectCloud1 = rectif(rectCloud1, 1,'Michelson',1) ;


rectCloud2 = random_cloud(rectZ2,1,0,0) ;
rectifiedRectCloud2 = rectif(rectCloud2, 1,'Michelson',1) ;


rectCloud3 = random_cloud(rectZ3,1,0,0) ;
rectifiedRectCloud3 = rectif(rectCloud3, 1,'Michelson',1) ;


rectCloud4 = random_cloud(rectZ4,1,0,0) ;
rectifiedRectCloud4 = rectif(rectCloud4, 1,'Michelson',1) ;


%%
figure
s1 = sum(z1,3) ;
s2 = sum(z2,3) ;
s3 = sum(z3,3) ;
s4 = sum(z4,3) ;
subplot(2,7,1)
plot(s1(:,ceil(N_X/2)))
hold on
plot(s2(:,ceil(N_X/2)),'r')
plot(s3(:,ceil(N_X/2)),'g')
plot(s4(:,ceil(N_X/2)),'k')

subplot(2,7,2)
imshow(rectifiedCloud1(:,:,1))
title(['frame with speed ',v1])

subplot(2,7,3)
imshow(rectifiedCloud2(:,:,2))
title(['frame with speed ',v2])

subplot(2,7,4)
imshow(rectifiedCloud3(:,:,3))
title(['frame with speed ',v3])
subplot(2,7,5)
imshow(rectifiedCloud4(:,:,4))
title(['frame with speed ',v4])

subplot(2,7,6)
d2 = s1./s2 ;
d3 = s1./s3 ;
d4 = s1./s4 ;
plot(d2(:,ceil(N_X/2)),'r')
hold on
plot(d3(:,ceil(N_X/2)),'g')
plot(d4(:,ceil(N_X/2)),'k')
title('ratio of spatial spectrum to v1');


subplot(2,7,7)
plot(s1(ceil(N_X/2),:))
hold on
plot(s2(ceil(N_X/2),:),'r')
plot(s3(ceil(N_X/2),:),'g')
plot(s4(ceil(N_X/2),:),'k')
title('spatial spectra');


rectZ1 = z1 ;
rectZ2 = z2./repmat(sum(z2,3),[1 1 , size(z1,3)]).*repmat(sum(z1,3),[1 1 , size(z1,3)]);
rectZ3 = z3./repmat(sum(z3,3),[1 1 , size(z1,3)]).*repmat(sum(z1,3),[1 1 , size(z1,3)]);
rectZ4 = z4./repmat(sum(z4,3),[1 1 , size(z1,3)]).*repmat(sum(z1,3),[1 1 , size(z1,3)]);

s1 = sum(rectZ1,3) ;
s2 = sum(rectZ2,3) ;
s3 = sum(rectZ3,3) ;
s4 = sum(rectZ4,3) ;
subplot(2,7,8)
plot(s1(:,ceil(N_X/2)))
hold on
plot(s2(:,ceil(N_X/2)),'r')
plot(s3(:,ceil(N_X/2)),'g')
plot(s4(:,ceil(N_X/2)),'k')

subplot(2,7,9)
imshow(rectifiedRectCloud1(:,:,1))

subplot(2,7,10)
imshow(rectifiedRectCloud2(:,:,2))

subplot(2,7,11)
imshow(rectifiedRectCloud3(:,:,3))

subplot(2,7,12)
imshow(rectifiedRectCloud4(:,:,4))

subplot(2,7,13)
d2 = s1./s2 ;
d3 = s1./s3 ;
d4 = s1./s4 ;
plot(d2(:,ceil(N_X/2)),'r')
hold on
plot(d3(:,ceil(N_X/2)),'g')
plot(d4(:,ceil(N_X/2)),'k')

subplot(2,7,14)
plot(s1(ceil(N_X/2),:))
hold on
plot(s2(ceil(N_X/2),:),'r')
plot(s3(ceil(N_X/2),:),'g')
plot(s4(ceil(N_X/2),:),'k')



%%
%%
% figure
% subplot(2,3,1)
% isosurface(z1)
% title(['spectrum, N =',num2str(N_frame1)])
% xlabel('fz')
% ylabel('fz')
% zlabel('fz')
% xlim([-0.5 0.5])
% ylim([-0.5 0.5])
% zlim([-0.5 0.5])
%%
% subplot(2,3,2)
% isosurface(env_rad1)
% isosurface(env_speed1)
% title('envelopes')
% xlabel('fz')
% ylabel('fz')
% zlabel('fz')
% xlim([-0.5 0.5])
% ylim([-0.5 0.5])
% zlim([-0.5 0.5])

% subplot(2,3,3)
% imshow(cloud1(:,:,1))
%
% subplot(2,3,4)
% isosurface(z1)
% title(['spectrum, N =',num2str(N_frame2)])
% xlabel('fz')
% ylabel('fz')
% zlabel('fz')
% xlim([-0.5 0.5])
% ylim([-0.5 0.5])
% zlim([-0.5 0.5])

% 
% subplot(2,3,5)
% isosurface(env_rad2)
% isosurface(env_speed2)
% title('envelopes')
% xlabel('fz')
% ylabel('fz')
% zlabel('fz')
% xlim([-0.5 0.5])
% ylim([-0.5 0.5])
% zlim([-0.5 0.5])




% subplot(2,3,6)
% imshow(cloud2(:,:,1))



