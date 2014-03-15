function vPixPerFrame = physicalSpeed2machineSpeed(vAnglePerSec, frameRate, viewingDistance, pixelPitch)
% from visual angle speed to pixel speed

pixPerAngle             = 2* viewingDistance * tand(1/2)/pixelPitch;
vPixPerSecond           = vAnglePerSec*pixPerAngle;
vPixPerFrame            = vPixPerSecond/frameRate;




