  
dataPath = '/export/data/motionClouds/MOCO_lifetime_MatFiles_30Hz_8.6deg_sec_MoreFrm' ;
allFilenames = {  'MOCO_random_F033_disp17_01.mat',
                  'MOCO_life002_F033_disp17_02.mat',
                  'MOCO_life003_F033_disp17_03.mat',
                  'MOCO_life004_F033_disp17_04.mat',
                  'MOCO_life005_F033_disp17_05.mat',
                  'MOCO_life006_F033_disp17_06.mat',
                  'MOCO_life007_F033_disp17_07.mat',
                  'MOCO_life008_F033_disp17_08.mat'};

  
  %%
%write to avi
for fileIdx = 1:length(allFilenames)%1:1%
  
  filename = allFilenames{fileIdx};
  load(fullfile(dataPath,filename))

  images =squeeze(images);  
  write2aviByYuv(uint32(images*255), fullfile(dataPath,strrep(filename,'.mat','.avi')), 30) ;
end