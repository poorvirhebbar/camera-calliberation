function [cameraParams] = estimateSingleCameraParameters(imagePoints, boardSize, patchSize, imageSize)
% This function will estimate camera parameters (intrinsic, extrinsic) from
% checkerboard image points.

% Zhang's method consists of 5 parts
% 1. Estimate homography from checkerboard plane to screen space.
% 2. Calculate B matrix by solving Vb = 0.
% 3. Extract intrinsic parameters from B matrix.
% 4. Calculate extrinsic parameters from intrinsic parameters and homography.
% 5. Refine parameters using the maximum likelihood estimation.

% Function inputs:
% - 'imagePoints': positions of checkerboard points in a screen space.
% - 'boardSize': the number of horizontal, vertical patchs in the checkerboard.
% - 'patchSize': the size of the checkerboard patch in mm.
% - 'imageSize': the size of the checkerboard image in pixels.

% Function outputs:
% - 'cameraParams': a camera parameter includes intrinsic and extrinsic.

numView = size(imagePoints, 3);
numVerticalPatch = boardSize(1) - 1;
numHorizontalPatch = boardSize(2) - 1;
numCorner = size(imagePoints, 1);

%% Estimate a homography (appendix A)
% Generate checkerboard world points
worldPoints = zeros(size(imagePoints,1), size(imagePoints,2));
for i=1:numHorizontalPatch
    for j=1:numVerticalPatch
        worldPoints(numVerticalPatch*(i-1)+j,1)=(i-1)*patchSize;
        worldPoints(numVerticalPatch*(i-1)+j,2)=(j-1)*patchSize;
    end
end

% Build L matrix
L = zeros(2 * numCorner, 9, numView);
X = zeros(numCorner,1);
Y = zeros(numCorner,1);
U = zeros(numCorner,numView);
V2= zeros(numCorner,numView);

for j=1:numCorner
    X(j)= worldPoints(j,1);
    Y(j)= worldPoints(j,2);
    for i=1:numView
        U(j,i)= imagePoints(j,1,i);
        V2(j,i)= imagePoints(j,2,i);
    end    
end

% Fill L matrix
for i = 1:numView
    for j = 1:numCorner  
        L(((2*j)-1),1,i) = -X(j);
        L(((2*j)-1),2,i) = -Y(j);
        L(((2*j)-1),3,i) = -1;
        L(((2*j)-1),7,i) = X(j)*U(j,i);
        L(((2*j)-1),8,i) = Y(j)*U(j,i);
        L(((2*j)-1),9,i) = U(j,i);
        L(2*j,4,i) = -X(j);
        L(2*j,5,i) = -Y(j);
        L(2*j,6,i) = -1;
        L(2*j,7,i) = X(j)*V2(j,i);
        L(2*j,8,i) = Y(j)*V2(j,i);
        L(2*j,9,i) = V2(j,i);
    end
end

homography = zeros(3,3,numView);
h = zeros(3,3,numView);
% Fill homography matrix

for i=1:numView
    [U1,S,V1] = svd(L(:,:,i)); 
    H=V1(:,9);
    H=H/H(9);
    homography(:,:,i)=transpose(reshape(H,3,3));
end
h = homography;
%% Solve closed-form (section 3.1)
h1= h(:,1,:);
h2= h(:,2,:);
h3= h(:,3,:);
%% Solve closed-form (section 3.1)
V = zeros(2 * numView, 6);
b = zeros(6, 1);


% Fill V matrix and calculate b vector
% ----- Your code here (4) ----- (slide 19, 23)
for i=1:numView
    v11 = [h(1,1,i)*h(1,1,i), (h(1,1,i)*h(2,1,i) + h(2,1,i)*h(1,1,i)), (h(1,1,i)*h(3,1,i) + h(3,1,i)*h(1,1,i)), h(2,1,i)*h(2,1,i), (h(2,1,i)*h(3,1,i) + h(3,1,i)*h(2,1,i)), h(3,1,i)*h(3,1,i)];
    v22 = [h(1,2,i)*h(1,2,i), (h(1,2,i)*h(2,2,i) + h(2,2,i)*h(1,2,i)), (h(1,2,i)*h(3,2,i) + h(3,2,i)*h(1,2,i)), h(2,2,i)*h(2,2,i), (h(2,2,i)*h(3,2,i) + h(3,2,i)*h(2,2,i)), h(3,2,i)*h(3,2,i)];
    v12 = [h(1,1,i)*h(1,2,i), (h(1,1,i)*h(2,2,i) + h(2,1,i)*h(1,2,i)), (h(1,1,i)*h(3,2,i) + h(3,1,i)*h(1,2,i)), h(2,1,i)*h(2,2,i), (h(2,1,i)*h(3,2,i) + h(3,1,i)*h(2,2,i)), h(3,1,i)*h(3,2,i)];
    V((2*i - 1),:) = v12;
    V(2*i,:)= v11 - v22;
end
[U1,S1,V1]=svd(V);
b = V1(:,6);

%% Extraction of the intrinsic parameters from matrix B (appendix B)

% ----- Your code here (5) ----- (slide 24)
B=zeros(3,3);
B(1,1)=b(1);
B(1,2)=b(2);
B(1,3)=b(3);
B(2,1)=b(2);
B(2,2)=b(4);
B(2,3)=b(5);
B(3,1)=b(3);
B(3,2)=b(5);
B(3,3)=b(6);

v0 = ((B(1,2)*B(1,3))-(B(1,1)*B(2,3)))/((B(1,1)*B(2,2))-(B(1,2)*B(1,2)));  

lambda = B(3,3) - ((((B(1,2)*B(1,3))-(B(1,1)*B(2,3)))*v0 + B(1,3)*B(1,3))/B(1,1));

alpha = sqrt(lambda/B(1,1)); 

beta = sqrt(lambda*B(1,1)/((B(1,1)*B(2,2))-(B(1,2)*B(1,2)))); 

gamma = -(B(1,2)*alpha*alpha*beta)/lambda;  

u0 = ((gamma*v0)/beta) - ((B(1,3)*(alpha*alpha))/lambda);  

p=[v0,lambda,alpha,beta,gamma,u0];
%% Estimate initial RT (section 3.1)
K = zeros(3,3);
K(1,1)=alpha;
K(1,2)=gamma;
K(1,3)=u0;
K(2,2)=beta;
K(2,3)=v0;
K(3,3)=1;
K1= inv(K);
z=zeros(numView,1);

for i=1:numView
    z=((1/(norm(K1*h1(:,i)))) + (1/(norm(K1*h2(:,i)))))/2;
    r1 = z*(K1*h1(:,i));
    r2 = z*(K1*h2(:,i));
    t  = z*(K1*h3(:,i));
    r3 = cross(r1,r2); 
    R = [r1,r2,r3];    
    [A,B,C] = svd(R);
    R1= A*(transpose(C));
    Rt(:,1:3,i)=R1;
    Rt(:,4,i)=t;
end
%% Maximum likelihood estimation (section 3.2)
options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', ...
    'TolX', 1e-32, 'TolFun', 1e-32, 'MaxFunEvals', 1e64, ...
    'MaxIter', 1e64, 'UseParallel', true);
x0 = zeros(5 + 6 * size(imagePoints, 3), 1); 
x0(1:5) = [alpha, beta, gamma, u0, v0];
rvec = zeros(3,numView);

for i=1:numView
    rvec(:,i)=rotationMatrixToVector(transpose(Rt(:,1:3,i)));    
end

for i=1:numView
    x0(6*i)    = Rt(1,4,i);
    x0(6*i + 1)= Rt(2,4,i);
    x0(6*i + 2)= Rt(3,4,i);
    x0(6*i + 3)= rvec(1,i);
    x0(6*i + 4)= rvec(2,i);
    x0(6*i + 5)= rvec(3,i);
end


% Non-least square optimization
% Read [https://mathworks.com/help/optim/ug/lsqnonlin.html] for more information
[objective] = @(x) func_calibration(imagePoints, worldPoints, x);

[x_hat, ~, ~, ~, ~] = lsqnonlin(objective,x0,[],[],options);


%% Build camera parameters
rvecs = zeros(numView, 3);
tvecs = zeros(numView, 3);
K = [1, 0, 0
     0, 1, 0
     0, 0, 1];

% Extract intrinsic matrix K, rotation vectors and translation vectors from x_hat
% ----- Your code here (8) -----
K(1,1)=x_hat(1);
K(2,2)=x_hat(2);
K(1,2)=x_hat(3);
K(1,3)=x_hat(4);
K(2,3)=x_hat(5);

for i=1:numView
    tvecs(i,1)=x_hat(6*i);
    tvecs(i,2)=x_hat(6*i + 1);
    tvecs(i,3)=x_hat(6*i + 2);
    rvecs(i,1)=x_hat(6*i + 3);
    rvecs(i,2)=x_hat(6*i + 4);
    rvecs(i,3)=x_hat(6*i + 5);    
end

% Generate cameraParameters structure
cameraParams = cameraParameters('IntrinsicMatrix', K', ...
    'RotationVectors', rvecs, 'TranslationVectors', tvecs, ...
    'WorldPoints', worldPoints, 'WorldUnits', 'mm', ...
    'imageSize', imageSize) ; 


reprojected_errors = zeros(size(imagePoints));

% Uncomment this line after you implement this function to calculate
% reprojection errors of your camera parameters.
 reprojected_errors = imagePoints - cameraParams.ReprojectedPoints;

cameraParams = cameraParameters('IntrinsicMatrix', K', ...
    'RotationVectors', rvecs, 'TranslationVectors', tvecs, ...
    'WorldPoints', worldPoints, 'WorldUnits', 'mm', ...
    'imageSize', imageSize, 'ReprojectionErrors', reprojected_errors) ; 