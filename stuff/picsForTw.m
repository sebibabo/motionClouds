addpath('functions')
addpath(genpath('/home/bosse/dev/matlab/0.1-bosse/functionPool'))
clear allv1              = 0;%physicalSpeed2machineSpeed(vAnglePerSec,frameRefresh,viewingDistance,pixelPitch) ;
v2              = 2;%physicalSpeed2machineSpeed(vAnglePerSec,frameRefresh,viewingDistance,pixelPitch) ;
v3              = 4;%physicalSpeed2machineSpeed(vAnglePerSec,frameRefresh,viewingDistance,pixelPitch) ;
v4              = 8;%physicalSpeed2machineSpeed(vAnglePerSec,frameRefresh,viewingDistance,pixelPitch) ;


%% set the the motion cloud parameters:
% stimulus dependent 


% hardware dependent:
N_X = 512; 
N_Y = N_X ;
N_frame = N_X/4 ;

alpha = 0 ;
theta = 0. ;
B_theta = inf ;
loggabor = true ;
ft_0 = 1;

[fx, fy, ft] = get_grids(N_X, N_Y, N_frame)  ;

sf_0  = 0.11 ;
B_sf  = 0.01;
Bv    = 0.11 ;
V_X = 4;
V_Y = 0;

B_V = 0.1;
env_speed = envelope_speed(fx, fy, ft, V_X, V_Y, B_V);
env_rad = envelope_radial(fx, fy, ft, sf_0, B_sf, loggabor) ;
%%    
isosurface(fx,fy,ft,env_speed);
isosurface(fx,fy,ft,env_rad);
xlabel('fx');
ylabel('fy');
zlabel('ft');
xlim([-0.5 0.5])
ylim([-0.5 0.5])
zlim([-0.5 0.5])

figure
env_speed = envelope_speed_adap(fx, fy, ft, V_X, V_Y, B_V);
env_rad = envelope_radial_adap(fx, fy, ft, sf_0, B_sf, loggabor) ;
    
isosurface(fx,fy,ft,env_speed);
isosurface(fx,fy,ft,env_rad);
xlabel('fx');
ylabel('fy');
zlabel('ft');
xlim([-0.5 0.5])
ylim([-0.5 0.5])
zlim([-0.5 0.5])

