clear all
close all

addpath(genpath('/home/bosse/dev/matlab/0.1-bosse/functionPool'))

directory = '/export/data/motionClouds/MOCO_coh_phaseRand_PowerDiva' ;

listing   = dir(fullfile(directory, 'MOCO_incohPhase_incohPhase*'));
offset    = 15 ;
imgIdx    = 3  ;
for fileIdx = 1:length(listing)
  
 scrsz = get(0,'ScreenSize');
%   
  
    figure('Position',[0 0  scrsz(3)/2 scrsz(4)/2]) ;
%   h = figure;
%  set(h,'PaperOrientation','landscape');
%   set(h, 'PaperUnits','normalized');
%   set(h, 'Position', [0 0 4 4]);



  load(fullfile(directory,listing(fileIdx).name));
  
  img_boil   = double(images(:,:,1,imgIdx));
  img_motion = double(images(:,:,1,imgIdx+offset));
  
  ampSpec_boil    = abs(fftshift(fft2(img_boil))) ;
  ampSpec_motion  = abs(fftshift(fft2(img_motion))) ;
  
  subplot(3,4,1)
  imshow(img_boil)
  axis image
  title(['Frame (phIncoh) ', num2str(imgIdx)]);
  subplot(3,4,2)
  imshow(img_motion)
  axis image
  title(['Frame (phIncoh) ', num2str(imgIdx+offset)]);
   
  subplot(3,4,5)
  imagesc(log(ampSpec_boil), [0, max(max([log(ampSpec_boil(:)),log(ampSpec_motion(:))]))])
  axis image
  title({'log-Spectrum',['frame (phIncoh) ', num2str(imgIdx)]});
  
  subplot(3,4,6)
  imagesc(log(ampSpec_motion),[0, max(max([log(ampSpec_boil(:)),log(ampSpec_motion(:))]))])
  axis image
  title({'log-Spectrum',['frame (phIncoh) ', num2str(imgIdx+offset)]});
  
  subplot(3,4,9)
  plot(sum(ampSpec_boil,1),'b');
  hold on
  plot(sum(ampSpec_boil,2),'r');
  legend('horizontal','vertical')
  title({'Spectrum along main axes',['frame (phIncoh) ', num2str(imgIdx)]});
  
  
  subplot(3,4,10)
  plot(sum(ampSpec_motion,1),'b');
  hold on
  plot(sum(ampSpec_motion,2),'r');

  title({'Spectrum along main axes',['frame (phIncoh) ', num2str(imgIdx+offset)]});
  %
  load(fullfile(directory,strrep(listing(fileIdx).name,'incohPhase_','')));
  
  img_boil   = double(images(:,:,1,imgIdx));
  img_motion = double(images(:,:,1,imgIdx+offset));
  
  ampSpec_boil    = abs(fftshift(fft2(img_boil))) ;
  ampSpec_motion  = abs(fftshift(fft2(img_motion))) ;
  
  subplot(3,4,3)
  imshow(img_boil)
  axis image
  title(['Frame (phCoh) ', num2str(imgIdx)]);
  subplot(3,4,4)
  imshow(img_motion)
  axis image
  title(['Frame (phCoh) ', num2str(imgIdx+offset)]);
  
  subplot(3,4,7)
  imagesc(log(ampSpec_boil), [0, max(max([log(ampSpec_boil(:)),log(ampSpec_motion(:))]))])
  axis image
  title({'log-Spectrum', ['frame (phCoh) ', num2str(imgIdx)]});
  
  subplot(3,4,8)
  imagesc(log(ampSpec_motion),[0, max(max([log(ampSpec_boil(:)),log(ampSpec_motion(:))]))])
  axis image
  title({'log-Spectrum', ['frame (phCoh) ', num2str(imgIdx+offset)]});
  
  subplot(3,4,11)
  plot(sum(ampSpec_boil,1),'b');
  hold on
  plot(sum(ampSpec_boil,2),'r');
  title({'Spectrum along main axes', ['frame (phCoh) ', num2str(imgIdx)]});
  
  
  subplot(3,4,12)
  plot(sum(ampSpec_motion,1),'b');
  hold on
  plot(sum(ampSpec_motion,2),'r');
  title({'Spectrum along main axes', ['frame (phCoh)', num2str(imgIdx+offset)]});
  
  suptitle(['Spectra of ',strrep(strrep(strrep(listing(fileIdx).name,'incohPhase_',''),'_','\_'),'.mat','')]);
  %saveTightFigure(h, fullfile(directory, ['cmpSpectra_',strrep(strrep(listing(fileIdx).name,'incohPhase_',''),'.mat',''),'.pdf']))
  
%   saveas(h,fullfile(directory,'temp.pdf'))
  export_fig( fullfile(directory, ['cmpSpectra_',strrep(strrep(listing(fileIdx).name,'incohPhase_',''),'.mat',''),'.pdf']))
  %%
end
close all

