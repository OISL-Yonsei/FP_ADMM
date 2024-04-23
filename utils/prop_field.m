function [out, TF] = prop_field(in,dx,lambda,propDist,RI,propMode)
% variables (all SI units):
% in - source plane field (2D or 1D)
% dx - pixel size (assuming square pixels for 2D)
% lambda - wavelength
% z - propagation distance
% mode 0: angular spectrum, mode 1: fresnel

if(propMode==0)
    [m,n] = size(in); 
    [fx,fy] = meshgrid((-m/2:m/2-1)./(m*dx), (-m/2:m/2-1)./(m*dx));
    k = 2*pi/lambda*RI;
    sqrtTerm=1-((fx.^2+fy.^2)*lambda.^2);
    sqrtTerm=sqrtTerm.*(sqrtTerm>0);
    TF = exp(1i*k*propDist*sqrt(sqrtTerm));
    ftIn = fftshift(fft2(in));
    out = ifft2(ifftshift(TF.*ftIn)); 
    sampling_criteria_fourier=(dx>=sqrt(lambda*propDist/m));
    
    if(sampling_criteria_fourier==0)
        fprintf('error(Nyquist rate), N should be more than %1f\n',lambda*propDist/(dx^2));
    end
    
elseif(propMode==1)
    % Input array size
    [m,n]=size(in); 
    %% calculation in fourier domain
    [fx,fy] = meshgrid((-m/2:m/2-1)./(m*dx), (-m/2:m/2-1)./(m*dx));
    k = 2*pi/lambda*RI;
    TF = exp(1i*k*propDist.*(1-(fx.^2+fy.^2).*lambda.^2./2));
    
    
    %% calculation in spatial domain
%     [x,y] = meshgrid((-m/2:m/2-1)./(dx), (-m/2:m/2-1)./(dx));
%     cpsf=exp((x.^2+y.^2).*1i.*k./(2.*propDist)).*exp(1i.*k.*propDist)./(1i.*lambda.*propDist);
%     TF=fftshift(fft2(cpsf));
    
    
    ftIn=fftshift(fft2(in));
    out = ifft2(ifftshift(TF.*ftIn)); 
    % Fourier transform of the convolution to the observation plane
    %Make sure aliasing
    sampling_criteria_fourier=(dx>=sqrt(lambda*propDist/m));
%     sampling_criteria_spatial=(dx<=(lambda*propDist/(m*dx));
    if(sampling_criteria_fourier==0)
        fprintf('error');
    end
    
end
