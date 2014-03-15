


function env_color = envelope_color(fx, fy, ft, alpha, ft_0)
%     Returns the color envelope. 
%     Run 'test_color.py' to see the effect of alpha
%     alpha = 0 white
%     alpha = 1 pink
%     alpha = 2 red/brownian
%     (see http//en.wikipedia.org/wiki/1/f_noise )
    ft_0=1;
    f_radius = frequency_radius(fx, fy, ft, ft_0).^alpha ;
    env_color = 1./f_radius ;
