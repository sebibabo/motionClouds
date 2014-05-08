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

v              = 4;%physicalSpeed2machineSpeed(vAnglePerSec,frameRefresh,viewingDistance,pixelPitch) ;
framesPerMod   = frameRefresh/stimModRate/2 ; 


%% set the the motion cloud parameters:
% stimulus dependent 

B_sf  = 0.01;
% hardware dependent:
N_X = 512; 
N_Y = N_X ;
N_frame1 = 512;% 2^size_T ;
N_frame2 = 512/4;% 2^size_T ;
N_frame3 = 512/8;% 2^size_T ;
N_frame4 = 512/16;% 2^size_T ;
V_X = v ;
V_Y = 0. ;

alpha = 0 ;
theta = 0. ;
B_theta = inf ;
  loggabor = true ;
ft_0 = 1;

if ceil(framesPerMod) ~= framesPerMod
    disp('frames per modulation is not integer')
end

[fx1, fy1, ft1] = get_grids(N_X, N_Y, N_frame1)  ;
[fx2, fy2, ft2] = get_grids(N_X, N_Y, N_frame2)  ;
[fx3, fy3, ft3] = get_grids(N_X, N_Y, N_frame3)  ;
[fx4, fy4, ft4] = get_grids(N_X, N_Y, N_frame4)  ;

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
sf_0 = 0.011 ;
Bv = 0.211 ;


%   name = 'rdk_Vx_'+num2str(V_X)+'_Bv1_'+ num2str(Bv_1)+ '_Bv2_' + num2str(Bv)+ '_Bsf_' + num2str(Bsf);
env_rad1 = envelope_radial_adap(fx1, fy1, ft1, sf_0, B_sf, loggabor) ;
env_speed1 =  envelope_speed_adap(fx1, fy1, ft1, V_X, V_Y, Bv) ;
z1 =  envelope_gabor_adap(fx1, fy1, ft1, V_X, V_Y,Bv, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
cloud1 = random_cloud(z1,1,0,0);
rectifiedZ1 = rectif(cloud1, 1,'Michelson',1) ;


env_rad2 = envelope_radial_adap(fx2, fy2, ft2, sf_0, B_sf, loggabor) ;
env_speed2 =  envelope_speed_adap(fx2, fy2, ft2, V_X, V_Y, Bv) ;
z2 =  envelope_gabor_adap(fx2, fy2, ft2, V_X, V_Y,Bv, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
cloud2 = random_cloud(z2,1,0,0) ;
rectifiedZ2 = rectif(cloud2, 1,'Michelson',1) ;

env_rad3 = envelope_radial_adap(fx3, fy3, ft3, sf_0, B_sf, loggabor) ;
env_speed3 =  envelope_speed_adap(fx3, fy3, ft3, V_X, V_Y, Bv) ;
z3 =  envelope_gabor_adap(fx3, fy3, ft3, V_X, V_Y,Bv, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
cloud3 = random_cloud(z3,1,0,0) ;
rectifiedZ3 = rectif(cloud3, 1,'Michelson',1) ;

env_rad4 = envelope_radial_adap(fx4, fy4, ft4, sf_0, B_sf, loggabor) ;
env_speed4 =  envelope_speed_adap(fx4, fy4, ft4, V_X, V_Y, Bv) ;
z4 =  envelope_gabor_adap(fx4, fy4, ft4, V_X, V_Y,Bv, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
cloud4 = random_cloud(z4,1,0,0) ;
rectifiedZ4 = rectif(cloud4, 1,'Michelson',1) ;

%%
figure
s1 = sum(z1,3) ;
s2 = sum(z2,3) ;
s3 = sum(z3,3) ;
s4 = sum(z4,3) ;
subplot(1,7,1)
plot(s1(:,ceil(N_X/2)))
hold on
plot(s2(:,ceil(N_X/2)),'r')
plot(s3(:,ceil(N_X/2)),'g')
plot(s4(:,ceil(N_X/2)),'k')

subplot(1,7,2)
imshow(rectifiedZ1(:,:,1))

subplot(1,7,3)
imshow(rectifiedZ2(:,:,2))

subplot(1,7,4)
imshow(rectifiedZ3(:,:,3))

subplot(1,7,5)
imshow(rectifiedZ4(:,:,4))

subplot(1,7,6)
d2 = s1./s2 ;
d3 = s1./s3 ;
d4 = s1./s4 ;
plot(d2(:,ceil(N_X/2)),'r')
hold on
plot(d3(:,ceil(N_X/2)),'g')
plot(d4(:,ceil(N_X/2)),'k')

subplot(1,7,7)
plot(s1(ceil(N_X/2),:))
hold on
plot(s2(ceil(N_X/2),:),'r')
plot(s3(ceil(N_X/2),:),'g')
plot(s4(ceil(N_X/2),:),'k')


