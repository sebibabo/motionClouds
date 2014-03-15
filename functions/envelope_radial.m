
function env_radial = envelope_radial(fx, fy, ft, sf_0, B_sf,  loggabor) % ft_0 is not used, but set to one
%     Radial frequency envelope
%     selects a sphere around a preferred frequency with a shell width B_sf.
%     Run 'test_radial.py' to see the explore the effect of sf_0 and B_sf
    if sf_0 == 0
      env_radial = 1 ;
    elseif loggabor
        [N_X, N_Y, N_frame] = size(fx) ;
        % see http//en.wikipedia.org/wiki/Log-normal_distribution
        fr = frequency_radius(fx, fy, ft, 1.) ;
        env_radial = 1./fr.*exp(-.5*(log(fr/sf_0).^2)/(log((sf_0+B_sf)/sf_0).^2)) ;
    else
        env_radial = exp(-.5*(frequency_radius(fx, fy, ft, 1.) - sf_0)^2/B_sf^2);
    end
