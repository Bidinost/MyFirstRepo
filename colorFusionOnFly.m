%% Environement configuration
%Load rgb and infrared camera's calibration parameters
load('stereoCalibParams.mat');
%Add vl-feat path
addpath('vlfeat-0.9.20');

%% Configure acquisition for color and depth
kinectHandles = mxNiCreateContextC320D320();
laserAcquisition = @mxNiDepthAcquisition;
colorAcquisition = @mxNiColorAcquisition;

f1=figure('KeyPressFcn',@spacePress);
srcDepthFrame=laserAcquisition(KinectHandles); 
srcDepthFrame = permute(srcDepthFrame,[2 1]);
srcColorFrame = mxNiColorAcquisition(KinectHandles); 
srcColorFrame=permute(srcColorFrame,[3 2 1]);
subplot(1,2,1),h1=imshow(srcColorFrame); 
subplot(1,2,2),h2=imshow(srcDepthFrame,[]);

colormap('jet');

spacePressed=0;

%% Realtime trajectory estimator
while (spacePressed~=1)
    tgtDepthFrame=laserAcquisition(KinectHandles); tgtDepthFrame=permute(tgtDepthFrame,[2 1]); 
    set(h2,'CDATA',tgtDepthFrame);

    tgtColorFrame=mxNiColorAcquisition(KinectHandles); tgtColorFrame=permute(tgtColorFrame,[3 2 1]);
    set(h1,'CDATA',tgtColorFrame);
    
    [R T] = colorFusionRototranslationEstimation(srcColorFrame,srcDepthFrame,tgtColorFrame,tgtDepthFrame);
    
    drawnow;
end

mxNiDeleteContext(kinectHandles);
close(f1);