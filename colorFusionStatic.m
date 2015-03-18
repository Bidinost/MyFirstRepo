%% Environement configuration
%Load rgb and infrared camera's calibration parameters
load('stereoCalibParams.mat');
%Load acquisitions
load('colorAndDepthAcquisition');
%Add vl-feat path
addpath('vlfeat-0.9.20');

%% Configure visualizer for color and depth

f1=figure('KeyPressFcn',@spacePress);
srcDepthFrame=depthMatrixes(:,:,1); 
srcColorFrame = rgbMatrixes(:,:,:,1); 
subplot(1,2,1),h1=imshow(srcColorFrame); 
subplot(1,2,2),h2=imshow(srcDepthFrame,[]);


%% Static trajectory estimator
for(i=2:size(depthMatrixes,3))
    tgtDepthFrame=depthMatrixes(:,:,i); 
    set(h2,'CDATA',tgtDepthFrame);

    tgtColorFrame=rgbMatrixes(:,:,:,i);
    set(h1,'CDATA',tgtColorFrame);
    
    [R T] = colorFusionRototranslationEstimation(srcColorFrame,srcDepthFrame,tgtColorFrame,tgtDepthFrame);
    
    drawnow;
end

close(f1);