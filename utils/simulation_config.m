% object parameter
resizePixel=600;
objPhaseScaleFactor=2*pi*1/4;
dpix_cam = 0.72e-6; % sampling pixel size of the CCD
dpix_obj = dpix_cam / 3; % final pixel size of the reconstruction


% System parameters
waveLength = 0.5e-6;
NA = 0.1;
radius=6;
LEDgap = 4e-3; % 4mm between adjacent LEDs
LEDheight = 50e-3; % distance bewteen the LED matrix and the sample
NA_syn=0.5;
nDiv=11.5;
z = 30e-6;

% Reconstruction parameters
Opts.maxIter=50;
Opts.step_max_obj = 1;
Opts.step_max_pupil= 0;
Opts.update_obj='local';% global/local
Opts.update_pupil='local';% global/local
Opts.loss_mode='amplitude';
Opts.step_mode='monotone'; % monotone, constant
Opts.adaptiveStep=0;
Opts.init_mode=1; % 1: interpolated, 2: random, 3: GT


%% Regularization Parameters
Opts.noise_mode='GaussianSNR_sameSTD'; %sameSTD
Opts.reg_obj=1;
Opts.reg_pupil=0;
Opts.nSubIter_obj=10;
Opts.nSubIter_pupil=10;
Opts.weight_regular_obj=0.01;     %0.005
Opts.weight_regular_pupil=5e-3; % 0.5
Opts.weight_split_obj=0.1;%1e-4;


%% ETC
Opts.adaptiveLowerBound=1;
Opts.adaptiveUpperBound=1;
Opts.vignettingFactor=1;
Opts.display=1;
Opts.mode='real';
Opts.loss_cal_mode='GS';
Opts.SNR=5;
SNR=Opts.SNR;