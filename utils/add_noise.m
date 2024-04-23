
function iMeasurements=add_noise(iMeasurements,Opts)
    nLed=size(iMeasurements);
    nLed=nLed(3);
    % maxVal_intensity_center=max(max(iMeasurements(:,:,(nLed+1)/2)));
    for iLED =1:nLed
        maxVal_intensity=max((max(iMeasurements(:,:,iLED))));
        varI = std2(iMeasurements(:,:,iLED)./maxVal_intensity)^2;
        noise_mean=0;
        if contains(Opts.noise_mode,'Gaussian')
            if contains(Opts.noise_mode,'sameSNR')
                noise_std = sqrt(varI/10^(Opts.SNR/10));
            elseif contains(Opts.noise_mode,'sameSTD')
                bf_img=iMeasurements(:,:,(nLed+1)/2);
                varI = std2(bf_img./(max(max(bf_img))))^2;
                noise_std = sqrt(varI/10^(Opts.SNR/10));
            end
            iMeasurements(:,:,iLED)=maxVal_intensity.*imnoise(iMeasurements(:,:,iLED)./maxVal_intensity,'gaussian',noise_mean,noise_std.^2);
            
            iMeasurements(:,:,iLED)=maxVal_intensity.*imnoise(iMeasurements(:,:,iLED)./maxVal_intensity,'poisson');
            
    
        elseif contains(Opts.noise_mode,'Poisson')
            iMeasurements(:,:,iLED)=maxVal_intensity.*imnoise(iMeasurements(:,:,iLED)./maxVal_intensity,'poisson');
        end
    
    end
end