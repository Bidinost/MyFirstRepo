%% Estimate rotation matrix and translation vector using SIFT for point matching
function [R, T] = colorFusionRototranslationEstimation(srcColorFrame, srcDepthFrame, tgtColorFrame, tgtDepthFrame, approach)
    [color2DMatches, scores] = getSiftMatches(srcColorFrame, tgtColorFrame);
    
    color2DMatches = getConsistentMatches(color2DMatches, scores, size(srcColorFrame));
    
    depth2DMatches = mapColorToIrCamera(color2DMatches, srcDepthFrame, tgtDepthFrame);
	
	[srcIrU, srcIrV] = ind2sub(size(srcDepthFrame),depth2DMatches(1,:));
	[tgtIrU, tgtIrV] = ind2sub(size(tgtDepthFrame),depth2DMatches(2,:));
	
	srcIrPoints = get3DPoints([srcIrU;srcIrV], srcDepthFrame);
	tgtIrPoints = get3DPoints([tgtIrU; tgtIrV], tgtDepthFrame);
	
	[R, T] = RTestimation(srcIrPoints, tgtIrPoints, approach, ones(3,size(srcIrPoints,2)) );
end

%Returns linearized indexes of hotspots in srcFrame and their match in
%tgtFrame
function [matches, scores] =  getSiftMatches(srcFrame, tgtFrame)
    srcFrameGray = single(rgb2gray(srcFrame));
    tgtFrameGray = single(rgb2gray(tgtFramme));
    
    %"frame" is a technical word in SIFT to define an hotspot
    [framesSrc, descriptorsSrc] = vl_sift(srcFrameGray) ;
    [framesTgt, descriptorsTgt] = vl_sift(tgtFrameGray) ;
    [matches, scores] = vl_ubcmatch(descriptorsSrc, descriptorsTgt);
    
    hotspotSrcU = round(framesSrc(1,matches(1,:)));
    hotspotSrcV = round(framesSrc(2,matches(1,:)));
    hotspotTgtU = round(framesTgt(1,matches(2,:)));
    hotspotTgtV = round(framesTgt(2,matches(2,:)));
    
    linearizedHotspotSrc = sub2ind(size(srcFrame),hotspotSrcU,hotspotSrcV);
    linearizedHotspotTgt = sub2ind(size(tgtFrame),hotspotTgtU,hotspotTgtV);
    
    matches = [linearizedHotspotSrc; linearizedHotspotTgt];
end

%Returns only "better" matches based on scores and pixel distance.
function [filteredMatches] = getConsistentMatches(matches, scores, frameSize)
    filteredMatches = [];
    pixelDistanceTresh = 10; % maximum permitted distance 
    scoreTresh = 50000; % maximum permitted descriptor's distance
    
    for (hotspotPair = [matches;scores])
        if(hotspotPair(3)>scoreTresh)
            continue;
        end
        
        [srcU, srcV] = ind2sub(frameSize, hotspotPair(1));
        [tgtU, tgtV] = ind2sub(frameSize, hotspotPair(2));
        
        if(norm([srcU;srcV] - [tgtU;tgtV]) <= pixelDistanceTresh)
            filteredMatches = [filteredMatches, [hotspotPair(1);hotspotPair(2)]];
        end
    end
end

function [irMatches] = mapColorToIrCamera(color2DMatches, srcIrFrame, tgtIrFrame)
	irSrcPointsLinearized = zeros(1,size(color2DMatches,2));
	irTgtPointsLinearized = zeros(1,size(color2DMatches,2));
	
	for(i=1:size(color2DMatches,2))
		[rgbSrcU, rgbSrcV] = ind2sub(size(srcIrFrame),color2DMatches(1,i));
		[rgbTgtU, rgbTgtV] = ind2sub(size(srcIrFrame),color2DMatches(2,i));
		
		[irSrcU, irSrcV] = convertToIr([rgbSrcU; rgbSrcV], srcIrFrame);
		[irTgtU, irTgtV] = convertToIr([rgbTgtU; rgbTgtV], tgtIrFrame);
		
		irSrcPointsLinearized(1,i) = sub2ind(size(srcIrFrame), irSrcU, irSrcV);
		irTgtPointsLinearized(2,i) = sub2ind(size(tgtIrFrame), irTgtU, irTgtV);
	end
	
	irMatches = [irSrcPointsLinearized; irTgtPointsLinearized];	
end

function [irU, irV] = convertToIr(rgbPixel, irImage)
	windowSearchSize = [21 21];
	[minU, maxU, minV, maxV] = searchBounds(rgbPixel, size(irImage), windowSearchSize);
	
	irCandidatePixels = [];
	distances = [];
	
	for(u=minU:maxU)
		for(v=minV:maxV)
			irCandidatePixels = [irCandidatePixels, [u,v]];
			
			%Calibration uses [v;u] instead of [u;v]
			z = single(irImage(u,v));
			irXYZ = [v*z; u*z; z];
			rgbXYZ = calib_KK_right*(calib_R*inv(calib_KK_left)*irXYZ + calib_T);
			rgbVU = round(rgbXYZ./rgbXYZ(3));
			rgbUV = [rgbVU(2);rgbVU(1)];
			
			distances = [distances, norm(rgbPixel - rgbUV)];
		end
	end
	
	[~, minIdx] = min(distances);
	[irU, irV] = irCandidatePixels(:,minIdx);
end

function [minU, maxU, minV, maxV] = searchBounds(rgbPixel, imageSize, windowSearchSize)
	halfHeigh = (windowSearchSize(1)-1)/2;
	halfWidth = (windowSearchSize(2)-1)/2;
	
	minU = max(1,rgbPixel(1)-halfHeigh);
	maxU = min(imageSize(1), rgbPixel(1)+halfHeigh);
	minV = max(1,rgbPixel(2)-halfWidth);
	maxV = max(imageSize(2),rgbPixel(2)+halfWidth);
end

function XYZ = get3DPoints(UV, depthImage)
	zU = UV(1).*UV(3);
	zV = UV(2).*UV(3);
	XYZ= inv(KK_left)*[zU; zV; UV(3)];
end