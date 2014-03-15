
  
  function env_gabor =  envelope_gabor_adap(fx, fy, ft, V_X, V_Y,...
                    B_V, sf_0, B_sf, loggabor,...
                    theta, B_theta, alpha)
    
%     Returns the Motion Cloud kernel
    env_color = 1;%envelope_color_adap(fx, fy, ft, alpha,1) ; 
    env_ori = envelope_orientation(fx, fy, ft, theta, B_theta)  ;
    env_rad = envelope_radial_adap(fx, fy, ft, sf_0, B_sf, loggabor) ;
    env_speed =  envelope_speed_adap(fx, fy, ft, V_X, V_Y, B_V) ;

    env_gabor = env_color .* env_ori .* env_rad.*  env_speed ;
