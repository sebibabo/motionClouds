function env_speed = envelope_speed_adap(fx, fy, ft, V_X, V_Y, B_V)
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

    env_speed = exp(-.5*((ft+fx*V_X+fy*V_Y)).^2./(B_V).^2); % ft_0 is 1!

