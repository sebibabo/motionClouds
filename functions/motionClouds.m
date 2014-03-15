
%%%%%
% Adpted from the python file by sebastian bosse
% Main script for generating Motion Clouds
% 
% (c) Laurent Perrinet - INT/CNRS

% ^functionault Motion Clouds parameters^
% 
% Motion Clouds (keyword) parameters
% size   	-- power of two to functionine the frame size (N_X, N_Y)
% size_T 	-- power of two to functionine the number of frames (N_frame)
% N_X 	-- frame size horizontal dimension [px]
% N_Y 	-- frame size vertical dimension [px]
% N_frame -- number of frames [frames] (a full period in time frames)
% alpha 	-- exponent for the color envelope.
% sf_0	-- mean spatial frequency relative to the sampling frequency.
% ft_0 	-- spatiotemporal scaling factor. 
% B_sf 	-- spatial frequency bandwidth
% V_X 	-- horizontal speed component
% V_Y	-- vertical speed component
% B_V	-- speed bandwidth
% theta	-- mean orientation of the Gabor kernel
% B_theta -- orientation bandwidth
% loggabor -- (boolean) if True it uses a logi-Gabor kernel 
% 
% Display parameters
% 
% vext   	-- movie format. Stimulus can be saved as a 3D (x-y-t) multimedia file .avi movie, .mat array, .zip folder with a frame sequence. 	
% ext    	-- frame image format.
% T_movie -- movie duration [s].
% fps  	-- frame per seconds

% functionault values
% size = 7
% size_T = 7
% N_X = 2^size
% N_Y = N_X
% N_frame = 2^size_T
% ft_0 = N_X/float(N_frame)
% alpha = 1.0 
% sf_0 = 0.15 
% B_sf = 0.1
% V_X = 1.
% V_Y = 0.
% B_V = .2
% theta = 0.
% B_theta = pi/32.
% loggabor = True
% vext = '.avi'
% ext =  '.png'
% T_movie = 8. % this value functionines the duration of a temporal period
% fps = int( N_frame / T_movie)


function [fx, fy, ft] =get_grids(N_X, N_Y, N_frame)
%     Use that function to functionine a reference outline for envelopes in Fourier space.
%     In general, it is more efficient to functionine dimensions as powers of 2.
    [fx, fy, ft] = meshgrid(floor(-N_X/2):floor((N_X-1)/2 + 1), floor(-N_Y/2):floor((N_Y-1)/2 + 1),floor(-N_frame/2):floor((N_frame-1)/2 + 1)) ;
    
    fx           = fx*1./N_X  ;
    fy           = fy*1./N_Y ;
    ft           = ft*1./N_frame ;

function fRadius = frequency_radius(fx, fy, ft, ft_0)
%      Returns the frequency radius. To see the effect of the scaling factor run
%      'test_color.py'
    [N_X, N_Y, N_frame] = size(fx) ;
    R2 = fx^2 + fy^2 + (ft/ft_0)^2 ; 
    R2(floor(N_X/2) , floor(N_Y/2) , floor(N_frame/2) ) = inf  ;
    fRadius = sqrt(R2) ;

function env_color = envelope_color(fx, fy, ft, alpha, ft_0)
%     Returns the color envelope. 
%     Run 'test_color.py' to see the effect of alpha
%     alpha = 0 white
%     alpha = 1 pink
%     alpha = 2 red/brownian
%     (see http//en.wikipedia.org/wiki/1/f_noise )
    f_radius = frequency_radius(fx, fy, ft, ft_0)^alpha ;
    env_color = 1./f_radius ;

function env_radial = envelope_radial(fx, fy, ft, sf_0, B_sf, ft_0, loggabor) % ft_0 is not used, but set to one
%     Radial frequency envelope
%     selects a sphere around a preferred frequency with a shell width B_sf.
%     Run 'test_radial.py' to see the explore the effect of sf_0 and B_sf
    if sf_0 == 0
      env_radial = 1 ;
    elseif loggabor
        [N_X, N_Y, N_frame] = size(fx) ;
        % see http//en.wikipedia.org/wiki/Log-normal_distribution
        fr = frequency_radius(fx, fy, ft, 1.) ;
        env_radial = 1./fr*exp(-.5*(log(fr/sf_0)^2)/(log((sf_0+B_sf)/sf_0)^2)) ;
    else
        env_radial = exp(-.5*(frequency_radius(fx, fy, ft, 1.) - sf_0)^2/B_sf^2);
    end
        

function env_speed = envelope_speed(fx, fy, ft, V_X, V_Y, B_V)
%      Speed envelope
%      selects the plane corresponding to the speed (V_X, V_Y) with some thickness B_V

%     (V_X, V_Y) = (0,1) is downward and  (V_X, V_Y) = (1,0) is rightward in the movie.
%      A speed of V_X=1 corresponds to an average displacement of 1/N_X per frame.
%      To achieve one spatial period in one temporal period, you should scale by
%      V_scale = N_X/float(N_frame)
%      If N_X=N_Y=N_frame and V=1, then it is one spatial period in one temporal
%      period. it can be seen in the MC cube. functionine ft_0 = N_X/N_frame
% 
%     Run 'test_speed.py' to explore the speed parameters

    
    [N_X, N_Y, N_frame] = size(fx);
    env_speed = exp(-.5*((ft+fx*V_X+fy*V_Y))^2/(B_V*frequency_radius(fx, fy, ft, 1))^2); % ft_0 is 1!

function envelope_dir = envelope_orientation(fx, fy, ft, theta, B_theta)
    
%     Orientation envelope
%     selects one central orientation theta, B_theta the spread
%     We use a von-Mises distribution on the orientation.
% 
%     Run 'test_orientation.py' to see the effect of changing theta and B_theta.
    
    if not(isinf( B_theta))
        angle = arctan2(fy, fx) ;
        envelope_dir = exp(cos(angle-theta)/B_theta) ;
        % along with its symmetric (because the output signal is real)
        envelope_dir = envelope_dir +exp(cos(angle-theta-pi)/B_theta) ;
        % and now selecting blobs
    else % for large bandwidth returns a strictly flat envelope
        envelope_dir = 1. ;
    end
    
  function env_gabor =  envelope_gabor(fx, fy, ft, V_X, V_Y,...
                    B_V, sf_0, B_sf, loggabor,...
                    theta, B_theta, alpha)
    
%     Returns the Motion Cloud kernel

    
    color = envelope_color(fx, fy, ft, alpha)
    env_gabor = color.*envelope_orientation(fx, fy, ft, theta, B_theta).*...
           envelope_radial(fx, fy, ft, sf_0, B_sf, loggabor).*...
           envelope_speed(fx, fy, ft, V_X, V_Y, B_V);

  function fRadius= frequency_radius_adap(fx, fy, ft, ft_0) % cylindrical envelope
%      Returns the frequency radius. To see the effect of the scaling factor run
%      'test_color.py'

    [N_X, N_Y, N_frame] = size(fx) ;
    R2 = fx^2 + fy^2   ;% cf . Paul Schrater 00
    R2(floor(N_X/2) , floor(N_Y/2) ,:)= inf  ;
    fRadius = sqrt(R2) ;

 % cylindrical envelope
  function env_radial = envelope_radial_adap(fx, fy, ft, sf_0, B_sf, ft_0, loggabor) % ft_0 is not used, but set to 1. 
%     Radial frequency envelope
%     selects a sphere around a preferred frequency with a shell width B_sf.
%     Run 'test_radial.py' to see the explore the effect of sf_0 and B_sf
    if sf_0 == 0
      env_radial = 1 ;
    elseif loggabor
        [N_X, N_Y, N_frame] = size(fx) ;
        % see http//en.wikipedia.org/wiki/Log-normal_distribution
        fr = frequency_radius_adap(fx, fy, 0, 1.) ;
        env_radial = 1./fr*exp(-.5*(log(fr/sf_0)^2)/(log((sf_0+B_sf)/sf_0)^2)) ;
    else
        env_radial = exp(-.5*(frequency_radius(fx, fy, ft, 1.) - sf_0)^2/B_sf^2);
    end
  
  function env_gabor =  envelope_gabor_adap(fx, fy, ft, V_X, V_Y,...
                    B_V, sf_0, B_sf, loggabor,...
                    theta, B_theta, alpha)
    
%     Returns the Motion Cloud kernel
    color = envelope_color_adap(fx, fy, ft, alpha) ; 
    env_ori = envelope_orientation(fx, fy, ft, theta, B_theta)  ;
    env_rad = envelope_radial_adap(fx, fy, ft, sf_0, B_sf, loggabor) ;
    env_speed =  envelope_speed(fx, fy, ft, V_X, V_Y, B_V) ;

    env = color * env_ori * env_rad*  env_speed ;

function z = random_cloud(envelope, seed, impulse, do_amp)
%     Returns a Motion Cloud movie as a 3D matrix.
%     It first creates a random phase spectrum and then it computes the inverse FFT to obtain
%     the spatiotemporal stimulus.
% 
%     - use a specific seed to specify the RNG's seed,
%     - test the impulse response of the kernel by setting impulse to True
%     - test the effect of randomizing amplitudes too by setting do_amp to True
    [N_X, N_Y, N_frame] = size(envelope);
    amps = 1.
    if impulse
        phase = 0.
    else
        rng(seed) ;
        phase = 2 * pi * rand(N_X, N_Y, N_frame) ;
        if do_amp
        disp('randomizing amps') ;
            amps = randn(N_X, N_Y, N_frame) ;
            % see Galerne, B., Gousseau, Y. & Morel, J.-M. Random phase textures Theory and synthesis. IEEE Transactions in Image Processing (2010). URL http//www.biomedsearch.com/nih/Random-Phase-Textures-Theory-Synthesis/20550995.html. (basically, they conclude "Even though the two processes ADSN and RPN have different Fourier modulus distributions (see Section 4), they produce visually similar results when applied to natural images as shown by Fig. 11.")
        end
    end
    Fz = amps.* envelope.*exp(1j * phase);
    % centering the spectrum
    Fz = ifftshift(Fz) ;
    Fz(0, 0, 0) = 0. ;
    z = real(ifftn((Fz))) ;
    


%%%%%%%%%%%%%%%%%%%%%%%%%% Display Tools %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function visualize(fx, fy, ft, z, azimuth=290., elevation=45., 
%     thresholds=[0.94, .89, .75, .5, .25,.1], opacities=[.9, .8, .7, .5, .2, .2],
%     name=None, ext=ext, do_axis=True, do_grids=False, draw_projections=True,
%     colorbar=False, f_N=2., f_tN=2., figsize=figsize)
%     
%      Visualize the  Fourier spectrum 
% 
%     N_X, N_Y, N_frame = z.shape
%     mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=figsize)
%     mlab.clf()
%     
%     % Normalize the amplitude.
%     z /= z.max() 
%     print "z.max = " + str(z.max() )
%     % Create scalar field
%     src = mlab.pipeline.scalar_field(fx, fy, ft, z)
%     if draw_projections
%         src_x = mlab.pipeline.scalar_field(fx,fy,ft,tile(sum(z,axis=0),(N_X,1,1)))
%         src_y = mlab.pipeline.scalar_field(fx,fy,ft,tile(reshape(sum(z,axis=1),(N_X,1,N_frame)), (1,N_Y,1)))
%         src_z = mlab.pipeline.scalar_field(fx,fy,ft,tile(reshape(sum(z,axis=2),(N_X,N_Y,1)), (1,1,N_frame)))
% 
%         % Create projections
%         border = 0.47
%         scpx = mlab.pipeline.scalar_cut_plane(src_x, plane_orientation='x_axes', view_controls=False)
%         scpx.implicit_plane.plane.origin = [-border, 1/N_Y, 1/N_frame]
%         scpx.enable_contours = True
%         scpy=mlab.pipeline.scalar_cut_plane(src_y,plane_orientation='y_axes', view_controls=False)
%         scpy.implicit_plane.plane.origin = [1/N_X, border, 1/N_frame]
%         scpy.enable_contours = True
%         scpz=mlab.pipeline.scalar_cut_plane(src_z,plane_orientation='z_axes', view_controls=False)
%         scpz.implicit_plane.plane.origin = [1/N_X, 1/N_Y, -border]
%         scpz.enable_contours = True
% 
%     % Generate iso-surfaces at differnet energy levels
%     for threshold, opacity in zip(thresholds, opacities)
%         mlab.pipeline.iso_surface(src, contours=[z.max()-threshold*z.ptp(), ],
%                                   opacity=opacity)
%         mlab.outline(extent=[-1./2,1./2,-1./2,1./2,-1./2,1./2],)
% 
%     % Draw a sphere at the origin
%     x = array([0])
%     y = array([0])
%     z = array([0])
%     s = 0.01
%     mlab.points3d(x, y, z, extent=[-s, s, -s, s, -s, s], scale_factor=0.15)
% 
%     if colorbar mlab.colorbar(title='density', orientation='horizontal')
%     if do_axis
%         ax = mlab.axes(xlabel='fx', ylabel='fy', zlabel='ft',
%                        extent=[-1./2, 1./2, -1./2, 1./2, -1./2, 1./2])
%         ax.axes.set(font_factor=2.)
%     
%     try
%         mlab.view(azimuth=azimuth, elevation=elevation, distance='auto',
%                   focalpoint='auto')
%     except
%         print(" You should upgrade your mayavi version")
% 
%     if not(name is None)
%         mlab.savefig(name + ext, magnification=1, size=figsize)
%     else
%        mlab.show(stop=True)
% 
%     mlab.close(all=True)

% function cube(fx, fy, ft, im, azimuth=-45., elevation=130., roll=-180., name=None,
%          ext=ext, do_axis=True, show_label=True, colormap='gray',
%          vmin=0., vmax=1., figsize=figsize)
% 
%     
%     Visualize the stimulus as a cube
%     
%     
% 
%     N_X, N_Y, N_frame = im.shape
%     mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=figsize)
%     mlab.clf()
%     src = mlab.pipeline.scalar_field(fx*2., fy*2., ft*2., im)
% 
%     mlab.pipeline.image_plane_widget(src, plane_orientation='z_axes',
%                                      slice_index=0, colormap=colormap, vmin=vmin, vmax=vmax)
%     mlab.pipeline.image_plane_widget(src, plane_orientation='z_axes',
%                                      slice_index=N_frame,colormap=colormap,
%                                      vmin=vmin, vmax=vmax)
%     mlab.pipeline.image_plane_widget(src, plane_orientation='x_axes', slice_index=0,
%                                      colormap=colormap, vmin=vmin, vmax=vmax)
%     mlab.pipeline.image_plane_widget(src, plane_orientation='x_axes', slice_index=N_X,
%                                      colormap=colormap, vmin=vmin, vmax=vmax)
% 
%     mlab.pipeline.image_plane_widget(src, plane_orientation='y_axes', slice_index=0,
%                                      colormap=colormap, vmin=vmin, vmax=vmax)
%     mlab.pipeline.image_plane_widget(src, plane_orientation='y_axes', slice_index=N_Y,
%                                      colormap=colormap, vmin=vmin, vmax=vmax)
% 
%     if do_axis
%         ax = mlab.axes(xlabel='x', ylabel='y', zlabel='t',
%                        extent=[-1., 1., -1., 1., -1., 1.], ranges=[0., N_X, 0., N_Y, 0., N_frame],
%                        x_axis_visibility=True, y_axis_visibility=True,
%                        z_axis_visibility=True)
%         ax.axes.set(font_factor=2.)
% 
%         if not(show_label) ax.axes.set(label_format='')
% 
%     
%     try
%         mlab.view(azimuth=azimuth, elevation=elevation, distance='auto',
%                   focalpoint='auto')
%         mlab.roll(roll=roll)
%     except
%         print(" You should upgrade your mayavi version")
%     
%     if not(name is None)
%         %mlab.options.offscreen = True
%         mlab.savefig(name + ext, magnification=1, size=figsize)
%     else
%         mlab.show(stop=True)
%     
%     mlab.close(all=True)




function z = rectif(z, contrast, method, verbose)
     
%     Transforms an image (can be 1,2 or 3D) with normal histogram into
%     a 0.5 centered image of determined contrast
%     method is either 'Michelson' or 'Energy'

    
    % Phase randomization takes any image and turns it into Gaussian-distributed noise of the same power (or, equivalently, variance).
    % See Peter J. Bex J. Opt. Soc. Am. A/Vol. 19, No. 6/June 2002 Spatial frequency, phase, and the contrast of natural images

    % Final rectification
    if verbose
        disp('Before Rectification of the frames')
        disp(['Mean=', num2str(mean(z(:))), ', std=', num2str(std(z(:))),', Min=',num2str(min(z(:))), ', Max=',num2str(max(z(:))), ' Abs(Max)=', num2str((max(abs(z(:)))))])
    end
    
    z = z-mean(z(:)); % this should be true *on average* in MotionClouds

    if strcmp(method,'Michelson')
        z = (.5* z/max(abs(z(:)))* contrast + .5) ;
    else
        z = (.5* z/std(z(:))  * contrast + .5) ;
    end
    if verbose
        disp( 'After Rectification of the frames')
        disp( ['Mean=', num2str(mean(z(:))), ', std=', num2str(std(z(:))), ', Min=', num2str(min(z(:))), ', Max=', num2str(max(z(:)))]);
        disp( ['percentage pixels clipped=',  num2str(sum(abs(z(:))>1.)*100/size(z))])
    end
    

