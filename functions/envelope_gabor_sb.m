    
  function env_gabor =  envelope_gabor_sb(fx, fy, ft, V_X, V_Y,...
                    B_V, sf_0, B_sf, loggabor,...
                    theta, B_theta, alpha)
    
%     Returns the Motion Cloud kernel

    
    color = envelope_color(fx, fy, ft, 0);
    env_gabor = color.*envelope_orientation(fx, fy, ft, theta, B_theta).*...
           envelope_radial(fx, fy, ft, sf_0, B_sf, loggabor).*...
           envelope_speed(fx, fy, ft, V_X, V_Y, B_V);

