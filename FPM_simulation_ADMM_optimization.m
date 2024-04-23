
%% 4/24/2024 Editted by KC Lee

% MIT License
% 
% Copyright (c) 2024 OISL-Yonsei
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE. 





%% simulate the forward imaging process of Fourier ptychography
% simulate the high resoluiton complex object
clear; clc; close all;
FT = @(x) fftshift(fft2((x))); %Fourier Transform
IFT = @(x) (ifft2(ifftshift(x))); %Inverse Fourier Transform

addpath('./utils');
addpath('./samples');
addpath('./config');

simulation_config();


%% simulate the high resoluiton complex object
testImgAmp = im2double(imread('westconcordorthophoto.png'));
testImgPhase = im2double(imread('./Samples/barbara.png'));

objAmp = imresize(testImgAmp,[resizePixel,resizePixel]);
objAmp = objAmp + mean2(objAmp) ./ 5;
objAmp = objAmp ./ max(max(objAmp));

objPhase = imresize(testImgPhase,[resizePixel, resizePixel]);
objPhase = objPhase - min(objPhase(:));
objPhase = objPhase ./ (max(max(objPhase)));
objPhase = objPhase - 0.5;
objPhase = objPhase .* objPhaseScaleFactor;


%% Simulated Object Generation
objBase= objAmp.*exp(1i.*objPhase);


%% setup the parameters for the coherent imaging system
k0 = 2*pi / waveLength;
nPixel_obj = resizePixel; % image size of high resolution object
FoV = nPixel_obj * dpix_obj;
du = 1 / FoV;


%% generate the low-pass filtered images

cutoffFrequency_syn = NA_syn * k0;
kmax_obj = pi / dpix_obj;
kx_unit_pos = 0:2*cutoffFrequency_syn/nDiv:cutoffFrequency_syn;
kx_unit = [-flip(kx_unit_pos(2:end)),kx_unit_pos];
[kxm_unit, kym_unit] = meshgrid(kx_unit,kx_unit);
criteria=((kxm_unit.^2+kym_unit.^2)<(((NA_syn+NA/2)*k0)^2));
nLed = sum(criteria(:));
Opts.nLed=nLed;
[row,col] = find(criteria);
for i =1:nLed
    kx(i)=kxm_unit(row(i),col(i));
    ky(i)=kym_unit(row(i),col(i));
end
kx_relative=kx./k0;
ky_relative=ky./k0;


% calculate parameters for k-space
nPixel_cam = nPixel_obj/(dpix_cam/dpix_obj); % image size of the final output
iMeasurements = zeros(nPixel_cam, nPixel_cam, nLed); % the final low resolution image sequence
dkx = 2*pi/(dpix_obj*nPixel_obj);
dky = 2*pi/(dpix_obj*nPixel_obj);
cutoffFrequency = NA * k0;
kmax = pi/dpix_cam;



[kxm, kym] = meshgrid(-kmax:kmax/((nPixel_cam-1)/2):kmax,-kmax:kmax/((nPixel_cam-1)/2):kmax);
CTF = double((kxm.^2+kym.^2)<cutoffFrequency^2); % pupil function circ(kmax)
cpsf = fftshift(ifft2(ifftshift(CTF))); % coherent PSF
ipsf = (abs(cpsf)).^2; % incoherent PSF
OTF = abs(fftshift(fft2(ifftshift(ipsf)))); % incoherent transfer function 
OTF = OTF./max(max(OTF));
kzm = sqrt(k0^2-kxm.^2-kym.^2);
pupil_defocus = exp(1i.*z.*real(kzm)).*exp(-abs(z).*abs(imag(kzm)));% defocus aberration
pupilGT = pupil_defocus.*CTF;
pupilGT = abs(pupilGT).*exp(1i.*(angle(pupilGT)-mean2(angle(pupilGT))));


%% Simulate intensity measurement 
objBaseFT = FT(objBase);
solutionSpace_obj = zeros(nPixel_obj);
solutionSpace_overlap = zeros(nPixel_obj);
overlapSpace = zeros(nPixel_obj);

cropR=zeros(nLed,4);
for iLED =1:nLed
    kxc = round((nPixel_obj+1)/2+kx(1,iLED)/dkx);
    kyc = round((nPixel_obj+1)/2+ky(1,iLED)/dky);
    kyl = round(kyc-(nPixel_cam-1)/2);
    kyh = round(kyc+(nPixel_cam-1)/2);
    kxl = round(kxc-(nPixel_cam-1)/2);
    kxh = round(kxc+(nPixel_cam-1)/2);
    cropR(iLED,1) = kyl; cropR(iLED,2) = kyh; cropR(iLED,3) = kxl; cropR(iLED,4) = kxh;
    lowResFT = objBaseFT(kyl:kyh,kxl:kxh).*pupilGT;
    iMeasurements(:,:,iLED) = (abs(IFT(lowResFT))).^2;
    solutionSpace_obj(kyl:kyh,kxl:kxh) = solutionSpace_obj(kyl:kyh,kxl:kxh)|CTF;
    solutionSpace_overlap(kyl:kyh,kxl:kxh) = solutionSpace_overlap(kyl:kyh,kxl:kxh) + CTF;
    if iLED==round((nLed+1)/2)|iLED==(round((nLed+1)/2)+1)
        overlapSpace(kyl:kyh,kxl:kxh)=overlapSpace(kyl:kyh,kxl:kxh)+CTF;
    end
end
r_overlap=sum(sum(overlapSpace==2))/sum((CTF(:)))*100;
fprintf('overlap is %.1f percent \n',r_overlap);



%% ADD noise to the measurement
iMeasurements = add_noise(iMeasurements,Opts);



%% Reorder for increasing convergence
illumNA=sqrt(kx_relative.^2+ky_relative.^2);
[dummy,seq]=sort(illumNA);
iMeasurements=iMeasurements(:,:,seq);
illumNA=illumNA(seq);
cropR=cropR(seq,:);

%% initialization

if Opts.init_mode==1 % interpolated initialization
    O0=mean2(iMeasurements);
    O0=imresize(O0,[nPixel_obj,nPixel_obj]);

elseif Opts.init_mode==2 % random initialization
    rng(1);
    O0=rand(nPixel_obj,nPixel_obj).*exp(1i.*rand(nPixel_obj,nPixel_obj));

elseif Opts.init_mode==3
    O0=objBase;
end

P0 = (pupilGT);
Ps = (CTF);
freqUV = [kx_relative,ky_relative];
Opts.nPixel_obj = nPixel_obj;
Opts.nPixel_cam = nPixel_cam;
Opts.kx = kx;
Opts.dkx = dkx;
Opts.ky = ky;
Opts.dky = dky;
Opts.O0 = O0;
Opts.P0 = P0;
Opts.Ps = Ps;
Opts.FT = FT;
Opts.IFT = IFT;
Opts.objGT = objBase;
Opts.pupilGT = pupilGT;
Opts.dfi = [];
Opts.H0 = [];
Opts.scale = [];
Opts.illumNA = sqrt(kx_relative.^2+ky_relative.^2);
Opts.NA = NA;
Opts.NA_obj = NA;
Opts.con = nPixel_cam.*dpix_obj;
Opts.cropR = cropR;
Opts.nIter = Opts.maxIter;
Opts.nImg = nLed;
Opts.r_overlap = r_overlap;


[highRes,pupil]=fpmFunc_optimization(iMeasurements,[nPixel_obj,nPixel_obj],freqUV,Opts);


