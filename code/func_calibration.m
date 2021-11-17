function [objective] = func_calibration(imagePoints, worldPoints, x)
% Objective function to minimize eq.10 in Zhang's paper. 
% Size of input variable x is 5+6*n where n is number of checkerboard 
% images. An intrinsic matrix can be reconstructed from first five
% parameters, and the extrinsic matrix can be reconstructed from remain
% parameters.

% You should fill the variable hat_m which contains reprojected positions 
% of checkerboard points in screen coordinate.

% Function inputs:
% - 'imagePoints': positions of checkerboard points in a screen space.
% - 'worldPoints': positions of checkerboard points in a model space.
% - 'x': parameters to be optimized.

% Function outputs:
% - 'objective': difference of estimated values and real values.
    
numView = size(imagePoints,3);
numCorner = size(imagePoints,1);
hat_m = zeros(size(imagePoints));

K = zeros(3,3);
K(1,1)=x(1);
K(2,2)=x(2);
K(1,2)=x(3);
K(1,3)=x(4);
K(2,3)=x(5);
K(3,3)=1;
for i=1:numView
    for j=1:numCorner
        R1 = rotationVectorToMatrix(x(6*i + 3 : 6*i + 5))';
        R=zeros(3,3);
        y=zeros(3);
        R(:,1:2)=R1(:,1:2);
        R(:,3)=x(6*i : 6*i + 2);
        y(1:2)=worldPoints(j,:);
        y(3)=1;
        z=(K*R)*y;
        hat_m(j,:,i)=z(1:2)/z(3);
        
    end
end
objective = imagePoints - hat_m;