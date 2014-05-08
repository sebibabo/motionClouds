addpath(genpath('/home/bosse/dev/matlab/0.1-bosse/functionPool'))
clear all
%% spectrum of moco
% dataPath = '/export/data/motionClouds/MOCO_lifetime_MatFiles';
% 
% allFilenames = {'MOCO_life002_rand_01.mat',
%                 'MOCO_life002_trans_02.mat',
%                 'MOCO_life003_trans_03.mat',
%                 'MOCO_life004_trans_04.mat',
%                 'MOCO_life005_trans_05.mat',
%                 'MOCO_life006_trans_06.mat',
%                 'MOCO_life007_trans_07.mat',
%                 'MOCO_life150_trans_08.mat',
%                 }
% 
% dataPath = '/export/data/motionClouds/MOCO_lifetime_MatFiles_moreFrm/';
% 
%  allFilenames = { 'MOCO_motion_F025_displ03_06.mat',
%                   'MOCO_motion_F025_displ12_04.mat',
%                   'MOCO_random_F025_displ03_05.mat',
%                     'MOCO_random_F025_displ12_03.mat'};
% 
% offset = 39 ;

dataPath = '/home/bosse/dev/motionClouds/mcForCmpWithMoco_coherence_highterSpatFreqOctBw_adap' ;
listing = dir(fullfile(dataPath,'*mat'));
for fileIdx = 1:length(listing)
  allFilenames{fileIdx} = listing(fileIdx).name;
end

                
                  
                  
                  
offset = 15 ;
for fileIdx = 1:length(allFilenames)%1:1%
  figure
  
  filename = allFilenames{fileIdx};
  load(fullfile(dataPath,filename))

  images =double(squeeze(images));

  part = 1;

  thisPartsImages = images(:,:,(part-1)*offset+1:(part-1)*offset+offset-1);
  fx = ([1:size(thisPartsImages,1)]-size(thisPartsImages,1)/2)/size(thisPartsImages,1);
  ft = ([1:size(thisPartsImages,3)]-size(thisPartsImages,3)/2)/size(thisPartsImages,3);
  
  
  spectrum = abs(fftshift(fftn(thisPartsImages)));

  subplot(3,4,1)
  temp = sum(spectrum,3);
  imagesc(fx,fx,(temp))
  title('spatial spectrum as sum across ft')
  xlabel('fx')
  ylabel('fy')
   
  
  subplot(3,4,2)
  plot(fx,(temp(floor(size(temp,1)/2),:)))
  title('intersection along fx = 0');
  ylabel('fy');
  
  subplot(3,4,5)
  temp = squeeze(sum(spectrum,2)) ;
  imagesc(fx,ft,(temp))
  title('spatiotemporal spectrum as sum across fy');
  xlabel('fx')
  ylabel('ft')
  
  subplot(3,4,6)
  plot(ft,(temp(floor(size(temp,1)/2),:)))
  title('intersection alog fx = 0');
  xlabel('ft')
  
  subplot(3,4,9)
  temp = squeeze(sum(spectrum,1)) ;
  imagesc(fx,ft,(temp))
  title('spatiotemporal spectrum as sum across fx');
  xlabel('fy')
  ylabel('ft')
    
  subplot(3,4,10)
  plot(ft,(temp(floor(size(temp,1)/2),:)))
  title('intersection alog fx = 0');
  xlabel('ft')
  %%
  part = 2;
  thisPartsImages = images(:,:,(part-1)*offset+1:(part-1)*offset+offset);
  fx = ([1:size(thisPartsImages,1)]-size(thisPartsImages,1)/2)/size(thisPartsImages,1);
  ft = ([1:size(thisPartsImages,3)]-size(thisPartsImages,3)/2)/size(thisPartsImages,3);
  
  spectrum = abs(fftshift(fftn(thisPartsImages)));

  subplot(3,4,3)
  temp = sum(spectrum,3);
  imagesc(fx,fx,(temp))
  title('spatial spectrum as sum across ft')
  xlabel('fx')
  ylabel('fy')
  
  subplot(3,4,4)
  plot(fx,(temp(floor(size(temp,1)/2),:)))
  title('intersection along fx = 0');
  ylabel('fy');

  subplot(3,4,7)
  temp = squeeze(sum(spectrum,2)) ;
  imagesc(fx,ft,(temp))
  title('spatiotemporal spectrum as sum across fy');
  xlabel('fx')
  ylabel('ft')
  
  subplot(3,4,8)
  plot(ft,(temp(floor(size(temp,1)/2),:)))
  title('intersection alog fx = 0');
  xlabel('ft')
  
  subplot(3,4,11)
  temp = squeeze(sum(spectrum,1)) ;
  imagesc(fx,ft,(temp))
  title('spatiotemporal spectrum as sum across fx');
  xlabel('fy')
  ylabel('ft')
   
  subplot(3,4,12)
  plot(ft,temp(ceil(floor(size(temp,1)/2)),:))
  title('intersection alog fx = 0');
  xlabel('ft')
  temporalSpectra{fileIdx} = temp(ceil(floor(size(temp,1)/2)),:);
  suplabel(strrep(filename,'_','\_'),'t');
end
%%
figure
temp = cell2mat(temporalSpectra')';
plot(ft,temp )
ylim([0, max(temp(:) )])
legend(strrep(allFilenames,'_','\_'))
title('temporal spectra as the interection at fx =0, after summing along fy') ;

