function z = random_cloud(envelope, seed, impulse, do_amp)
%     Returns a Motion Cloud movie as a 3D matrix.
%     It first creates a random phase spectrum and then it computes the inverse FFT to obtain
%     the spatiotemporal stimulus.
% 
%     - use a specific seed to specify the RNG's seed,
%     - test the impulse response of the kernel by setting impulse to True
%     - test the effect of randomizing amplitudes too by setting do_amp to True
    [N_X, N_Y, N_frame] = size(envelope);
    amps = 1.;
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
    Fz(1, 1, 1) = 0. ;
    z = real(ifftn((Fz))) ;
    