

function [fx, fy, ft] =get_grids(N_X, N_Y, N_frame)
%     Use that function to functionine a reference outline for envelopes in Fourier space.
%     In general, it is more efficient to functionine dimensions as powers of 2.
    [fx, fy, ft] = meshgrid(floor(-N_X/2):floor((N_X-1)/2), floor(-N_Y/2):floor((N_Y-1)/2 ),floor(-N_frame/2):floor((N_frame-1)/2)) ;
    
    fx           = fx*1./N_X  ;
    fy           = fy*1./N_Y ;
    ft           = ft*1./N_frame ;

