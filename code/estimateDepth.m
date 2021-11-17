function [depthMap, disparityMap] = estimateDepth(leftImage, rightImage, stereoParameters)
% This function estimate disparity and depth values from left and right
% images. You should calculate disparty map first and then convert the
% disparity map to depth map using left camera parameters.

% Function inputs:
% - 'leftImage': rectified left image.
% - 'rightImage': rectified right image.
% - 'stereoParameters': stereo camera parameters.

% Function outputs:
% - 'depth': depth map of left camera.
% - 'disparity': disparity map of left camera.

leftImageGray = rgb2gray(im2double(leftImage));
rightImageGray = rgb2gray(im2double(rightImage));

translation = stereoParameters.TranslationOfCamera2;
baseline = norm(translation);
focalLength = stereoParameters.CameraParameters1.FocalLength(1);

disparityMap = zeros(size(leftImageGray));
depthMap = zeros(size(leftImageGray));

max_d = 583;%assuming width to be maximum disparity
min_d = 0;
w_radius=5;
kernel = ones(w_radius*2+1);
kernel = kernel ./ numel(kernel);

d_vals = min_d : max_d;
num_d = length(d_vals);
C = NaN(size(leftImageGray,1), size(leftImageGray,2), num_d); % the cost volume

for i = 1 : length(d_vals)
    d = d_vals(i);
    I2t = imtranslate(rightImageGray, [d 0]);
    C(:,:,i) = abs(leftImageGray - I2t);
    C(:,:,i) = imfilter(C(:,:,i), kernel);    
end

[C_min, D] = min(C, [], 3);
D = D + min_d;
a=size(D,1);
b=size(D,2);
disparityMap=D;
for i=1:a
    for j=1:b
        depthMap(i,j) = focalLength*baseline/D(i,j);
    end
end
end

