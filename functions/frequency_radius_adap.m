
  function fRadius= frequency_radius_adap(fx, fy, ft, ft_0) % cylindrical envelope
%      Returns the frequency radius. To see the effect of the scaling factor run
%      'test_color.py'

    [N_X, N_Y, N_frame] = size(fx) ;
    R2 = fx.^2 + fy.^2   ;% cf . Paul Schrater 00
    R2(floor(N_X/2)+1 , floor(N_Y/2)+1 ,:)= inf  ;
    fRadius = sqrt(R2) ;
