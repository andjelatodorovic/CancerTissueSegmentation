function [ ProbabilityMaps ] = EstUsingGMM( Image )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
k = 2;

BgMask = BgEstimate(Image, 0.87);

luma = double(rgb2gray(Image));
stainsLuma = luma(BgMask==0);


%otsu
thresh = multithresh(stainsLuma,2);
lbl = imquantize(luma,thresh);

%iv = reshape(double(Image(BgMask)), [],3);
%iv = reshape(rgb2hsv(double(Image)),[],3);
%{
sinV = sin(iv(:,1)*2*pi);
cosV = cos(iv(:,1)*2*pi);
newIV = zeros(size(iv,1), 3);
%newIV(:,1) = sinV/2;
newIV(:,1) = cosV/2;
newIV(:,2) = iv(:,2);
newIV(:,3) = iv(:,3);
%}

%{
iv = applycform(Image, makecform('srgb2lab'));
iv = reshape(double(iv), [],3);
%}

%{
lambda = 0.05;
bestVal = 10000000;
numPix = size(iv,1);

for i = 3:5 
    GMModel = fitgmdist(iv,i,'Options',statset('Display','final','MaxIter',200,'TolFun',1e-6));
    [idxvgmm,nlogl,P] = cluster(GMModel,iv);
    val = nlogl / numPix + lambda*i;

    if (val < bestVal)
        bestVal = val;
        k = i;
    end    
end
%}

%GMModel = fitgmdist(stainsLuma,k,'Regularize',0, 'Options',statset('Display','final','MaxIter',1500,'TolFun',1e-6));
%load('gmmmodel.mat');
%[idxvgmm,nlogl,P] = cluster(GMModel,stainsLuma);
%P = reshape(P, size(Image,1),size(Image,2),k);

%yuv = rgb2ycbcr(GMModel.mu);

%{
mu = zeros(k,3);
angle = atan2(GMModel.mu(:,1), GMModel.mu(:,2))
mask = angle < 0; 
angle(mask) = angle(mask) + 2*pi;
mu(:,1) = angle / (2*pi);
mu(:,2) = GMModel.mu(:,3);
mu(:,3) = GMModel.mu(:,4);
yuv = rgb2ycbcr(hsv2rgb(mu));
%}

%yuv = rgb2ycbcr(applycform(GMModel.mu, makecform('lab2srgb')));
%Luma = yuv(:,1);

% Luma = GMModel.mu;
% [maxC, bgId] = max(Luma);
% [minC, nucleiId] = min(Luma);

bgProb = zeros(size(Image,1),size(Image,2));
bgProb(BgMask) = 1.0;
nucleiProb = zeros(size(Image,1),size(Image,2));
%nucleiProb(BgMask==0) = P(:,nucleiId);
nucleiProb(lbl == 1) = 1.0;
ProbabilityMaps = zeros(size(Image));
ProbabilityMaps(:,:,1) = bgProb;%P(:,:,bgId);
ProbabilityMaps(:,:,3) = nucleiProb;%P(:,:,nucleiId);
ProbabilityMaps(:,:,2) = 1.0 - ProbabilityMaps(:,:,1) - ProbabilityMaps(:,:,3);
§
end

