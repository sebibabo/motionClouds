
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
