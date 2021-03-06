clear all
close all
addpath(genpath('/home/bosse/dev/matlab/0.1-bosse/functionPool'))
addpath(genpath('/home/bosse/dev/motionClouds/functions'))
%% spectrum of moco
dataPath = '/export/data/motionClouds/MOCO_coh_MatFiles_30Hz_moreFrm';


% 
allFilenames = {'MOCO_coh000_F033_displ17_01.mat',
  'MOCO_coh010_F033_displ17_02.mat',
  'MOCO_coh020_F033_displ17_03.mat',
  'MOCO_coh050_F033_displ17_04.mat',
  'MOCO_coh100_F033_displ17_05.mat',
                  }

listing = dir(fullfile(dataPath,'*mat'));
% 
% for fileIdx = 1:length(listing)
%   allFilenames{fileIdx} = listing(fileIdx).name;
% end
% 
                
                  
                  
                  
offset = 45 ;
for fileIdx = 1:length(allFilenames)%1:1%
  figure
  
  filename = allFilenames{fileIdx};
  load(fullfile(dataPath,filename))

  images =squeeze(images);

  part = 1;

  thisPartsImages = images(:,:,(part-1)*offset+1:(part-1)*offset+offset);
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
  plot(fx,(temp(size(temp,1)/2,:)))
  title('intersection along fx = 0');
  ylabel('fy');
  
  subplot(3,4,5)
  temp = squeeze(sum(spectrum,2)) ;
  imagesc(fx,ft,(temp))
  title('spatiotemporal spectrum as sum across fy');
  xlabel('fx')
  ylabel('ft')
  
  subplot(3,4,6)
  plot(ft,(temp(size(temp,1)/2,:)))
  title('intersection alog fx = 0');
  xlabel('ft')
  
  subplot(3,4,9)
  temp = squeeze(sum(spectrum,1)) ;
  imagesc(fx,ft,(temp))
  title('spatiotemporal spectrum as sum across fx');
  xlabel('fy')
  ylabel('ft')
    
  subplot(3,4,10)
  plot(ft,(temp(size(temp,1)/2,:)))
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
  plot(fx,(temp(size(temp,1)/2,:)))
  title('intersection along fx = 0');
  ylabel('fy');

  subplot(3,4,7)
  temp = squeeze(sum(spectrum,2)) ;
  imagesc(fx,ft,(temp))
  title('spatiotemporal spectrum as sum across fy');
  xlabel('fx')
  ylabel('ft')
  
  subplot(3,4,8)
  plot(ft,(temp(size(temp,1)/2,:)))
  title('intersection alog fx = 0');
  xlabel('ft')
  
  subplot(3,4,11)
  temp = squeeze(sum(spectrum,1)) ;
  imagesc(fx,ft,(temp))
  title('spatiotemporal spectrum as sum across fx');
  xlabel('fy')
  ylabel('ft')
   
  subplot(3,4,12)
  plot(ft,temp(ceil(size(temp,1)/2),:))
  title('intersection alog fx = 0');
  xlabel('ft')
  temporalSpectra{fileIdx} = temp(ceil(size(temp,1)/2),:);
  suplabel(strrep(filename,'_','\_'),'t');
end
%%
figure

plot(ft,cell2mat(temporalSpectra')')
legend(strrep(allFilenames,'_','\_'))
title('temporal spectra as the interection at fx =0, after summing along fy') ;

%% get fwhm
for fileIdx = 1:length(allFilenames)
  fwhmFreq_temp{fileIdx}   = getFwhm(temporalSpectra{fileIdx});
  amplitude_ratio{fileIdx} = max(temporalSpectra{fileIdx})/min(temporalSpectra{fileIdx}) ;

end



for fileIdx = 1:length(allFilenames)
  
  spectrum = temporalSpectra{fileIdx} ;
  [peakVal,peakIdx] =max(spectrum );
  % left side first
  leftSpectrum = spectrum (1:peakIdx-1);
  [leftVal,leftIdx] = min(abs(leftSpectrum-peakVal/2));
  leftVal           = leftSpectrum(leftIdx);
  
  if leftVal < peakVal/2 % half value between leftIdx and leftIdx+1;
    thisLeftFt = ft(leftIdx) ;
    thisRightFt = ft(leftIdx+1) ;
    
    thisLeftValue = spectrum(leftIdx);
    thisRightValue = spectrum(leftIdx+1);
    
  elseif leftVal > peakVal/2 % half value between leftIdx-1 and leftIdx;
    if leftIdx ==1 % inf to the left
      thisLeftFt  = inf ;
      thisRightFt = ft(leftIdx) ;
      
      thisLeftValue  = ft(leftIdx)   ;
      thisRightValue = spectrum(leftIdx) ;
    else
      thisLeftFt = ft(leftIdx-1) ;
      thisRightFt = ft(leftIdx) ;
      
      thisLeftValue  = ft(leftIdx-1) ;
      thisRightValue = spectrum(leftIdx) ;
    end
  end
  
  % linear interpolation
  m = (thisRightValue-thisLeftValue)/(thisRightFt-thisLeftFt);
  b = (thisLeftValue*thisRightFt - thisRightValue*thisLeftFt)/(thisRightFt-thisLeftFt);
  
  leftFreq = 1/m*(peakVal/2-b);
  
  % now the right side
  rightSpectrum = spectrum (peakIdx+1:end);
  [rightVal,rightIdx] = min(abs(rightSpectrum-peakVal/2));
  rightVal            = rightSpectrum(rightIdx) ;
  % bring rightIdx to full spectrum:
  rightIdx            = rightIdx+peakIdx;
  if spectrum(rightIdx) ~= rightVal
    error('something fishy');
  end
  
  if rightVal < peakVal/2 % half value between rightIdx-1 and rightIdx;
    thisLeftFt = ft(rightIdx-1) ;
    thisRightFt = ft(rightIdx) ;
    
    thisLeftValue = spectrum(rightIdx-1);
    thisRightValue = spectrum(rightIdx);
    
    
  elseif rightVal> peakVal/2 % half value between rightIdx and leftIdx+1;
    if rightIdx  == length(spectrum) % inf to the right
      thisLeftFt  = inf ;
      thisRightFt = ft(rightIdx) ;
      
      thisLeftValue  = ft(rightIdx)   ;
      thisRightValue = spectrum(rightIdx) ;
    else
      thisLeftFt = ft(rightIdx) ;
      thisRightFt = ft(rightIdx+1) ;
      
      thisLeftValue  = ft(rightIdx) ;
      thisRightValue = spectrum(rightIdx+1) ;
    end
  end
  % linear interpolation
  m = (thisRightValue-thisLeftValue)/(thisRightFt-thisLeftFt);
  b = (thisLeftValue*thisRightFt - thisRightValue*thisLeftFt)/(thisRightFt-thisLeftFt);
  
  rightFreq = 1/m*(peakVal/2-b);
  fwhmFreq{fileIdx} = rightFreq-leftFreq;
end

