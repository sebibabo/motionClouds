




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
        disp( ['percentage pixels clipped=',  num2str(sum(abs(z(:))>1.)*100/numel(z))])
    end
    z(z>1) =1 ;
    z(z<0) =0 ;

