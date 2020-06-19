frame178 = imread('frame178.tif');
%frame178 = rgb2gray(frame178);
frame178 = double(frame178);
frame179 = imread('frame179.tif');
frame179 = double(frame179);


motionEstimation(frame178, frame179, 30);
