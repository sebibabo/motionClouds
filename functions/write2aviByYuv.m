function write2aviByYuv(matrix, filename, freq)
%%
  
[m n t] = size(matrix);
  dummyU  = 127*ones(m/2,n/2,t);
  
  dummyU  = mat2cell(dummyU, m/2, n/2, ones(1,t));
   
  B=mat2cell(matrix, m, n, ones(1,t));
  yuv_export(B, dummyU, dummyU, 'temp.yuv', length(B), 'w');

  command = ['ffmpeg -y -s ', num2str(size(matrix,1)), 'x', num2str(size(matrix,1)), ' -r ', num2str(freq), ' -i temp.yuv  -vcodec rawvideo ' filename];
  system(command);
  command = 'rm temp.yuv';
  system(command);