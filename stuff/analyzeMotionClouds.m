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
N_frame = N_X/4;


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
sf_0 = 0.015;
B_sf = 0.001 ;
  
v  = 0 ;
Bv = 0.5 ;
paramList{1} = [v, Bv];

v  = 8 ;
Bv = 0.0001 ;
paramList{2} = [v, Bv];


for paramsIdx = 1:length(paramList)
  v = paramList{paramsIdx}(1);
  Bv = paramList{paramsIdx}(2);
  
  %   name = 'rdk_Vx_'+num2str(V_X)+'_Bv1_'+ num2str(Bv_1)+ '_Bv2_' + num2str(Bv)+ '_Bsf_' + num2str(Bsf);
  z{paramsIdx} =  envelope_gabor_adap(fx, fy, ft, v, V_Y,Bv, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
  cloud        = random_cloud(z{paramsIdx},1,0,0);
  rectifiedClouds{paramsIdx} = rectif(cloud, 1,'Michelson',1) ;
  s{paramsIdx} = sum(z{paramsIdx},3) ;  
end
%%
for paramsIdx = 1:length(paramList)
  subplot(3,2,paramsIdx)
  imshow(rectifiedClouds{paramsIdx}(:,:,1)) 
  title(['1st frame of cloud, [v, Bv] = ',num2str(paramList{paramsIdx})])
  
  subplot(3,2,paramsIdx+2)
  tempHist = rectifiedClouds{paramsIdx}(:,:,1);
  hist(tempHist(:),255)
end
%%
subplot(3,2,5)
plot(fx(1,:,1),s{1}(:,ceil(N_X/2)))
hold on
plot(fx(1,:,1),s{2}(:,ceil(N_X/2)),'r')
title('Projection to spatial frequencies of 1st frame')
legend({['[v, Bv] = ',num2str(paramList{1})],['[v, Bv] = ',num2str(paramList{2})]});

subplot(3,2,6)
plot(fx(1,:,1),s{1}(ceil(N_X/2),:))
hold on
plot(fx(1,:,1),s{2}(ceil(N_X/2),:),'r')

%%
for paramsIdx = 1:length(paramList)
  implay(rectifiedClouds{paramsIdx});
end
%%
figure
subplot(1,2,1)
plot(fx(1,:,1),s{1}(:,ceil(N_X/2))./s{2}(:,ceil(N_X/2)))
title('Ration between spatial frequencies ')

subplot(1,2,2)
plot(fx(1,:,1),s{1}(ceil(N_X/2),:)./s{2}(ceil(N_X/2),:))
title('Ration between spatial frequencies ')

%%
% isosurface(fx,fy,ft,z{1})
% xlim([-0.5 0.5])
% ylim([-0.5 0.5])
% zlim([-0.5 0.5])
% figure
% isosurface(fx,fy,ft,z{2})
% xlim([-0.5 0.5])
% ylim([-0.5 0.5])
% zlim([-0.5 0.5])

%%
% 1. possibility: dont care
% 2. possibility: incoherent: just randomize the frames
% 3. possibility: normalize spatial frequencies
  
%   name = 'rdk_Vx_'+num2str(V_X)+'_Bv1_'+ num2str(Bv_1)+ '_Bv2_' + num2str(Bv)+ '_Bsf_' + num2str(Bsf);
scaledZ{1}    = z{1};
scaledZ{2}    = z{2}./repmat(sum(z{2},3),[1,1,size(z{1},3)]).*repmat(sum(z{1},3),[1,1,size(z{1},3)]);
scaledZ{2}(isnan(scaledZ{2})) = 0;
cloud         = random_cloud(scaledZ{1},1,0,0);
rectifiedScaledClouds{1} = rectif(cloud, 1,'Michelson',1) ;
cloud         = random_cloud(scaledZ{2},1,0,0);
rectifiedScaledClouds{2} = rectif(cloud, 1,'Michelson',1) ;


scaledS{1} = sum(scaledZ{1},3) ;  
scaledS{2} = sum(scaledZ{2},3) ;  

%%
for paramsIdx = 1:length(paramList)
  subplot(3,2,paramsIdx)
  imshow(rectifiedScaledClouds{paramsIdx}(:,:,1)) 
  title(['1st frame of scaled cloud, [v, Bv] = ',num2str(paramList{paramsIdx})])
  subplot(3,2,paramsIdx+2)
  hist(rectifiedScaledClouds{paramsIdx}(:),255)
end
%%
subplot(3,2,5)
plot(fx(1,:,1),scaledS{1}(:,ceil(N_X/2)))
hold on
plot(fx(1,:,1),scaledS{2}(:,ceil(N_X/2)),'r')
title('Projection to spatial frequencies of 1st frame')
legend({['[v, Bv] = ',num2str(paramList{1})],['[v, Bv] = ',num2str(paramList{2})]});

subplot(3,2,6)
plot(fx(1,:,1),scaledS{1}(ceil(N_X/2),:))
hold on
plot(fx(1,:,1),scaledS{2}(ceil(N_X/2),:),'r--')

%%
for paramsIdx = 1:length(paramList)
  implay(rectifiedScaledClouds{paramsIdx});
end
%%
figure
subplot(1,2,1)
plot(fx(1,:,1),scaledS{1}(:,ceil(N_X/2))./scaledS{2}(:,ceil(N_X/2)))
title('Ration between spatial frequencies ')

subplot(1,2,2)
plot(fx(1,:,1),scaledS{1}(ceil(N_X/2),:)./scaledS{2}(ceil(N_X/2),:))
title('Ration between spatial frequencies ')

% %%
% isosurface(fx,fy,ft,scaledZ{1})
% xlim([-0.5 0.5])
% ylim([-0.5 0.5])
% zlim([-0.5 0.5])
% figure
% isosurface(fx,fy,ft,scaledZ{2})
% xlim([-0.5 0.5])
% ylim([-0.5 0.5])
% zlim([-0.5 0.5])

%% 4. histeq
for paramsIdx = 1:length(paramList)
   
  temp = reshape(rectifiedScaledClouds{paramsIdx}, size(rectifiedClouds{paramsIdx},1),[]) ;
  h    = length(temp(:))*1/255*ones(1,255);
  temp = histeq(temp);
  histeqClouds{paramsIdx} = reshape( temp, size(rectifiedClouds{paramsIdx},1),size(rectifiedClouds{paramsIdx},2),size(rectifiedClouds{paramsIdx},3));
  
  histeqZ{paramsIdx}      = abs(fftshift(fftn(histeqClouds{paramsIdx}) ));
  
  histeqS{paramsIdx} = sum(histeqZ{paramsIdx},3) ;  
end

%%
figure
for paramsIdx = 1:length(paramList)
  subplot(3,2,paramsIdx)
  imshow(histeqClouds{paramsIdx}(:,:,1)) 
  title(['1st frame of cloud, [v, Bv] = ',num2str(paramList{paramsIdx})])
  
  subplot(3,2,paramsIdx+2)
  hist(histeqClouds{paramsIdx}(:),255)
end
%%
subplot(3,2,5)
plot(fx(1,:,1),histeqS{1}(:,ceil(N_X/2)))
hold on
plot(fx(1,:,1),histeqS{2}(:,ceil(N_X/2)),'r')
title('Projection to spatial frequencies of 1st frame')
legend({['[v, Bv] = ',num2str(paramList{1})],['[v, Bv] = ',num2str(paramList{2})]});

subplot(3,2,6)
plot(fx(1,:,1),histeqS{1}(ceil(N_X/2),:))
hold on
plot(fx(1,:,1),histeqS{2}(ceil(N_X/2),:),'r--')

%%
for paramsIdx = 1:length(paramList)
  implay(histeqClouds{paramsIdx});
end
%%
figure
subplot(1,2,1)
plot(fx(1,:,1),histeqS{1}(:,ceil(N_X/2))./s{2}(:,ceil(N_X/2)))
title('Ration between spatial frequencies ')

subplot(1,2,2)
plot(fx(1,:,1),histeqS{1}(ceil(N_X/2),:)./s{2}(ceil(N_X/2),:))
title('Ration between spatial frequencies ')

%%
%  another solution: adapt energies
for paramsIdx = 1:length(paramList)
  v = paramList{paramsIdx}(1);
  Bv = paramList{paramsIdx}(2);
  
  %   name = 'rdk_Vx_'+num2str(V_X)+'_Bv1_'+ num2str(Bv_1)+ '_Bv2_' + num2str(Bv)+ '_Bsf_' + num2str(Bsf);
  energyScaledZ{paramsIdx} =  envelope_gabor_adap2(fx, fy, ft, v, V_Y,Bv, sf_0, B_sf, loggabor,theta, B_theta, alpha) ;
  cloud        = random_cloud(z{paramsIdx},1,0,0);
  energyScaledRectifiedClouds{paramsIdx} = rectif(cloud, 1,'Michelson',1) ;
  energyScaledS{paramsIdx} = sum(energyScaledZ{paramsIdx},3) ;  
end
%%
for paramsIdx = 1:length(paramList)
  subplot(3,2,paramsIdx)
  imshow(energyScaledRectifiedClouds{paramsIdx}(:,:,1)) 
  title(['1st frame of cloud, [v, Bv] = ',num2str(paramList{paramsIdx})])
  
  subplot(3,2,paramsIdx+2)
  hist(energyScaledRectifiedClouds{paramsIdx}(:),255)
end
%%
subplot(3,2,5)
plot(fx(1,:,1),energyScaledS{1}(:,ceil(N_X/2)))
hold on
plot(fx(1,:,1),energyScaledS{2}(:,ceil(N_X/2)),'r')
title('Projection to spatial frequencies of 1st frame')
legend({['[v, Bv] = ',num2str(paramList{1})],['[v, Bv] = ',num2str(paramList{2})]});

subplot(3,2,6)
plot(fx(1,:,1),energyScaledS{1}(ceil(N_X/2),:))
hold on
plot(fx(1,:,1),energyScaledS{2}(ceil(N_X/2),:),'r')

%%
for paramsIdx = 1:length(paramList)
  implay(energyScaledRectifiedClouds{paramsIdx});
end
%%
figure
subplot(1,2,1)
plot(fx(1,:,1),energyScaledS{1}(:,ceil(N_X/2))./energyScaledS{2}(:,ceil(N_X/2)))
title('Ration between spatial frequencies ')

subplot(1,2,2)
plot(fx(1,:,1),energyScaledS{1}(ceil(N_X/2),:)./energyScaledS{2}(ceil(N_X/2),:))
title('Ration between spatial frequencies ')
