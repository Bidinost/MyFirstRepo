%% Estimate rotation matrix and translation vector using SIFT for point matching
function [R, T] = colorFusionRototranslationEstimation(srcColorFrame, srcDepthFrame, tgtColorFrame, tgtDepthFrame, approach)
    [color2DMatches, scores] = getSiftMatches(srcColorFrame, tgtColorFrame);
    
    color2DMatches = getConsistentMatches(color2DMatches, scores, size(srcColorFrame));
    
    depth2DMatches = mapColorToIrCamera(color2DMatches, srcDepthFrame, tgtDepthFrame);
	
	if(isempty(depth2DMatches))
		R = nan;
		T = nan;
		return;
	end
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
	irMatches = [];
	
	for(colorPointPair=color2DMatches)
		[rgbSrcU, rgbSrcV] = ind2sub(size(srcIrFrame),colorPointPair(1));
		[rgbTgtU, rgbTgtV] = ind2sub(size(srcIrFrame),colorPointPair(2);
		
		[irSrcU, irSrcV] = convertToIr([rgbSrcU; rgbSrcV], srcIrFrame);
		[irTgtU, irTgtV] = convertToIr([rgbTgtU; rgbTgtV], tgtIrFrame);
		
		if(~isnan(irSrcU*irTgtU))
			irMatches = [irMatches, [sub2ind(size(srcIrFrame), irSrcU, irSrcV); sub2ind(size(tgtIrFrame), irTgtU, irTgtV)]];
		end
	end
end

function [irU, irV] = convertToIr(rgbPixel, irImage)
	windowSearchSize = [11 11];
	nearTresh = 5; %Maximum distance for valid mapping
	
	[minU, maxU, minV, maxV] = searchBounds(rgbPixel, size(irImage), windowSearchSize);
	
	Ugrid = repmat((minU:maxU)',1,maxV-minV+1);
	Vgrid = repmat((minV:maxV), maxU-minU+1,1);
	Zgrid = single(irImage(minV:maxV, minU:maxU));
	Zgrid(Zgrid==0)=nan;
	
	%%{
	% Without for loop
	irXYZ = [Vgrid(:)'.*Zgrid(:)'; Ugrid(:)'.*Zgrid(:)'; Zgrid(:)'];
	rgbXYZ = calib_KK_right*(calib_R*inv(calib_KK_left)*irXYZ + repmat(calib_T,1,size(irXYZ,2)));
	rgbVU = round(rgbXYZ./repmat(rgbXYZ(3),3,1));
	rgbUV = [rgbVU(2,:); rgbVU(1,:)]
	
	distances = sqrt(sum((rgbUV - repmat(rgbPixel,1,size(irXYZ,2))).^2,2));
	
	[minValue, minIdx] = min(distances);
	%If there's no valid match or rgb mapped point is too far from original rgb point, return nan
	if(isnan(minValue) || minValue > nearTresh)
		irU = nan;
		irV = nan;
		return;
	end
	
	irU = rgbUV(1,minIdx);
	irV = rgbUV(2,minIdx);
	%}
	
	%{
	% With for loop
	irCandidatePixels = [];
	distances = [];
	
	for(u=minU:maxU)
		for(v=minV:maxV)
			%Calibration uses [v;u] instead of [u;v]
			z = single(irImage(u,v));
			%Exclude saturated depth pixels
			if(z==0)
				continue;
			end
			
			irXYZ = [v*z; u*z; z];
			rgbXYZ = calib_KK_right*(calib_R*inv(calib_KK_left)*irXYZ + calib_T);
			rgbVU = round(rgbXYZ./rgbXYZ(3));
			rgbUV = [rgbVU(2);rgbVU(1)];
			
			irCandidatePixels = [irCandidatePixels, [u,v]];
			distances = [distances, norm(rgbPixel - rgbUV)];
		end
	end
	
	if(isempty(irCandidatePixels))
		irU = nan;
		irV = nan;
		return;
	end
	
	[minValue, minIdx] = min(distances);
	%If rgb mapped point is too far from original rgb point, return nan
	if(minValue > nearTresh)
		irU = nan;
		irV = nan;
		return;
	end
	
	irU = irCandidatePixels(1,minIdx);
	irV = irCandidatePixels(2,minIdx);
	%}
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