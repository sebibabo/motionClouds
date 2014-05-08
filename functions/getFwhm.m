function fwhm = getFwhm(spectrum)

  ft = [1:length(spectrum)]/length(spectrum)-0.5;
  [peakVal,peakIdx] =max(spectrum );
  % left side first
  leftSpectrum = spectrum (1:peakIdx-1);
  [leftVal,leftIdx] = min(abs(leftSpectrum-peakVal/2));
  leftVal           = leftSpectrum(leftIdx);
  
  if all(leftSpectrum>peakVal/2)
    leftFreq = -1 ;
  else
    if leftVal < peakVal/2 % half value between leftIdx and leftIdx+1;
      thisLeftFt = ft(leftIdx) ;
      thisRightFt = ft(leftIdx+1) ;

      thisLeftValue = spectrum(leftIdx);
      thisRightValue = spectrum(leftIdx+1);

    elseif leftVal > peakVal/2 % half value between leftIdx-1 and leftIdx;
      if leftIdx ==1 % inf to the left
        thisLeftFt  = 0.75 ;
        thisRightFt = ft(leftIdx) ;

        thisLeftValue  = ft(leftIdx)   ;
        thisRightValue = spectrum(leftIdx) ;
      else
        thisLeftFt = ft(leftIdx-1) ;
        thisRightFt = ft(leftIdx) ;

        thisLeftValue  = ft(leftIdx-1) ;
        thisRightValue = spectrum(leftIdx) ;
      end
    end

    % linear interpolation
    m = (thisRightValue-thisLeftValue)/(thisRightFt-thisLeftFt);
    b = (thisLeftValue*thisRightFt - thisRightValue*thisLeftFt)/(thisRightFt-thisLeftFt);

    leftFreq = 1/m*(peakVal/2-b);
  end
  % now the right side
  rightSpectrum = spectrum (peakIdx+1:end);
  [rightVal,rightIdx] = min(abs(rightSpectrum-peakVal/2));
  rightVal            = rightSpectrum(rightIdx) ;
  % bring rightIdx to full spectrum:
  rightIdx            = rightIdx+peakIdx;
  if spectrum(rightIdx) ~= rightVal
    error('something fishy');
  end
  
  if all(rightSpectrum>peakVal/2)
    rightFreq = 1;
  else
    if rightVal < peakVal/2 % half value between rightIdx-1 and rightIdx;
      thisLeftFt = ft(rightIdx-1) ;
      thisRightFt = ft(rightIdx) ;

      thisLeftValue = spectrum(rightIdx-1);
      thisRightValue = spectrum(rightIdx);


    elseif rightVal> peakVal/2 % half value between rightIdx and leftIdx+1;
      if rightIdx  == length(spectrum) % inf to the right
        thisLeftFt  = 0.75 ;
        thisRightFt = ft(rightIdx) ;

        thisLeftValue  = ft(rightIdx)   ;
        thisRightValue = spectrum(rightIdx) ;
      else
        thisLeftFt = ft(rightIdx) ;
        thisRightFt = ft(rightIdx+1) ;

        thisLeftValue  = ft(rightIdx) ;
        thisRightValue = spectrum(rightIdx+1) ;
      end
    end
    % linear interpolation
    m = (thisRightValue-thisLeftValue)/(thisRightFt-thisLeftFt);
    b = (thisLeftValue*thisRightFt - thisRightValue*thisLeftFt)/(thisRightFt-thisLeftFt);

    rightFreq = 1/m*(peakVal/2-b);
  end
  fwhm = rightFreq-leftFreq;

end