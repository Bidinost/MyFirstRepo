function [R, T] = RTestimation(srcPoints,targetPoints,approach, weights)

	switch approach
		case 'leastSquareSmall'
			[R, T] = leastSquareSmallEstimation(srcPoints, targetPoints);
		case 'leastSquareNonLinear'
			[R, T] = leastSquareNonLinearEstimation(srcPoints, targetPoints);
		case 'SVD'
			[R, T] = svdEstimation(srcPoints, targetPoints, weights);
	end
end

%Assuming small rototraslation between point clouds, solving Ax=b
%
% that is, denoting points as (x,y,z), 'source' as 's', 'target' as 't', 
% theta angles as 'th', traslation as 'T':
%
% xs - thz*ys + thy*zs + Tx = xt
% thz*xs + ys - thx*zs + Ty = yt
% -thy*xs + thx*ys + zs + Tz = zt
%
% that is:
%
% |  0    zs   -ys   1   0   0 |      | thx |      | xt |
% | -zs   0     xs   0   1   0 |   *  | thy |   =  | yt |
% |  ys  -xs    0    0   0   1 |      | thz |      | zt |
%                                     | Tx  |
%                                     | Ty  |
%                                     | Tz  |
%
% related to a generic source-target point pair
function [R, T] = leastSquareSmallEstimation(srcPoints, targetPoints)
	n = size(srcPoints,2);
	
	reshepableConstants = [zeros(1,n);         %0
							srcPoints(3,:);    %zs
							-srcPoints(2,:);   %-ys
							ones(1,n);         %1 
							zeros(2,n);        %0; 0
							-srcPoints(3,:);   %-zs
							zeros(1,n);        %0
							srcPoints(1,:);    %xs
							zeros(1,n);        %0
							ones(1,n);         %1
							zeros(1,n);        %0
							srcPoints(2,:);    %ys
							-srcPoints(1,:);   %-xs
							zeros(3,n);        %0; 0; 0;
							ones(1,n)];        %1
							
	A = reshape(reshepableConstants,3*n,6)';
	
	b = targetPoints(:);
	
	W = diag(weights);
	
	x = inv(A'*A)*A'*b;
	
	R = [1 -x(3) x(2);
		x(3) 1 -x(1);
		-x(2) x(1) 1];
	
	T = x(4:6);
end

%For a generic rototranslation between point cloud, solve Ax=b, starting
%from generic rotation and traslation:
%
% | cos(thx)*cos(thy) cos(thz)*sin(thy)*sin(thx)-sin(thz)*cos(thx)	cos(thz)*sin(thy)*cos(thx)+sin(thz)*sin(thx) |   | xs |   | Tx |   | xt |
% | sin(thz)*cos(thy) sin(thz)*sin(thy)*sin(thx)+cos(thz)*cos(thx)	sin(thz)*sin(thy)*cos(thx)-cos(thz)*sin(thx) | * | ys | + | Ty | = | yt |
% | -sin(thy)		  cos(thy)*sin(thx)								cos(thy)*cos(thx)						     |   | zs |   | Tz |   | zt |
%
% that is, denoting points as (x,y,z), 'source' as 's', 'target' as 't',
% traslation as 'T' and other unknow parameters with letters {A,..,I}:
%
% A*xs + B*ys + C*zs + Tx = xt
% D*xs + E*ys + F*zs + Ty = yt
% G*xs + H*ys + I*zs + Tz = zt
%
% that is:
%
% | xs ys zs 0  0  0  0  0  0  1  0  0 |      | A  |      | xt |
% | 0  0  0  xs ys zs 0  0  0  0  1  0 |   *  | B  |   =  | yt |
% | 0  0  0  0  0  0  xs ys zs 0  0  1 |      | C  |      | zt |
%                                             | D  |
%                                             | E  |
%                                             | F  |
%                                             | G  |
%                                             | H  |
%                                             | I  |
%                                             | Tx |
%                                             | Ty |
%                                             | Tz |
function [R, T] = leastSquareNonLinearEstimation(srcPoints, targetPoints)
	n = size(srcPoints,2);
	
	reshepableConstants= [srcPoints;           %xs; ys; zs;
							zeros(6,n);        %0; 0; 0; 0; 0; 0;
							ones(1,n);         %1
							zeros(5,n);        %0; 0; 0; 0; 0;
							srcPoints;         %xs; ys; zs;
							zeros(4,n);        %0; 0; 0; 0;
							ones(1,n);         %1
							zeros(7,n);        %0; 0; 0; 0; 0; 0; 0;
							srcPoints;         %xs; ys; zs;
							zeros(2,n);        %0; 0;
							ones(1,n)  ];      %1
							
	A = reshape(reshepableConstants,3*n,12)';
	
	b = targetPoints(:);
	
	W = diag(weights);
	
	x = inv(A'*A)*A'*b;
	
	R = [x(1) x(2) x(3);
		x(4) x(5) x(6);
		x(7) x(8) x(9) ];
	
	T = x(10:12);
end

function [R, T] = svdEstimation(srcPoints, targetPoints, weights)
	nSrcPoints = size(srcPoints,2);
	nTgtPoints = size(targetPoints,2);

	% normalize weights
	weights = weights ./ sum(weights);

	% find data centroid and deviations from centroid
	targetPoints_bar = targetPoints * transpose(weights);
	targetPoints_mark = targetPoints - repmat(targetPoints_bar, 1, nTgtPoints);
	% Apply weights
	targetPoints_mark = targetPoints_mark .* repmat(weights, 3, 1);
	
	% find data centroid and deviations from centroid
	srcPoints_bar = srcPoints * transpose(weights);
	srcPoints_mark = srcPoints - repmat(srcPoints_bar, 1, nSrcPoints);

	
	N = srcPoints_mark*transpose(targetPoints_mark);

	[U,~,V] = svd(N); % singular value decomposition

	R = V*diag([1 1 det(U*V')])*transpose(U);

	T = targetPoints_bar - R*srcPoints_bar;
end