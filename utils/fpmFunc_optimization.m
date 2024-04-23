function [highRes, pupil, loss, scale, freqUV] = fpmFunc_optimization( I, nObj, freqUV, Opts)

%% Default options
if nargin < 4
    disp('check reconstruction config');
end

%% Operators & derived constants
nSampR = size(I); % size of measurement
nImg = nSampR(3);
nSampR(3) = [];
cen0 = round((nObj+1)/2); % coordinate of center of Fourier space

%% Operators
row = @(x) x(:).';
% operator to crop region of O from proper location at the O plane
downsamp = @(x,cen) x(cen(2)-floor(nSampR(1)/2):cen(2)-floor(nSampR(1)/2)+nSampR(1)-1,...
    cen(1)-floor(nSampR(2)/2):cen(1)-floor(nSampR(2)/2)+nSampR(2)-1);
% cen is in (x,y)
poscost = @(ss) sum(sum((sqrt(abs(IFT(downsamp(highResFT, round(ss)).*pupil.*H0)).^2)-sqrt(I_mea)).^2));

%% Initialization
% copy some parameters
Ps = Opts.Ps;
H0 = Opts.H0;
FT = Opts.FT;
IFT = Opts.IFT;
dfi=Opts.dfi;
illumNA=Opts.illumNA;
NA=Opts.NA;
pupil = Opts.P0; 
%initialize pupil to GT
% pupil = Opts.pupilGT; 


Opts.P0 = 0; %Pupil function in Fourier domain
highRes = Opts.O0; Opts.O0 = 0; %Object in real space
highResFT=FT(highRes); %Init in Fourier space
con=Opts.con;
adaptiveLowerBound=Opts.adaptiveLowerBound;
adaptiveUpperBound=Opts.adaptiveUpperBound;
lossPrevious = inf;
lossCurrent = 50;
loss = zeros(1,Opts.maxIter);
iter = 0;
scale = Opts.scale;
step_max_obj = Opts.step_max_obj;
step_max_pupil= Opts.step_max_pupil;
step_pupil=zeros(1,Opts.maxIter);
step_obj=zeros(1,Opts.maxIter);
energy_obj=zeros(Opts.maxIter,4);


weight_regular_obj=Opts.weight_regular_obj;     %0.01
weight_regular_pupil=Opts.weight_regular_pupil; % 0.5
weight_split_obj=Opts.weight_split_obj;


%% Main
% convert from 1/um center to corresponding crop region
pupilShiftY = round(freqUV(:,2)*con);
pupilShiftX = round(freqUV(:,1)*con);

% corresponding crop regions
YXmid=floor(nObj./2)+1;
hfSz=floor(nSampR(1)/2);
cropR=Opts.cropR;
Ps=(Ps);
U      = (zeros(nObj));
U_0    = (zeros(nObj));
U_1    = (zeros(nObj));
V0     = (zeros(nObj));
V1     = (zeros(nObj));
W0     = (zeros(nObj));
W1     = (zeros(nObj));
I=(I);
pupil=(pupil);

% Generate pupil constraint for calculating laplacian
diamond_kernel=strel('diamond',1);
Ps_erode=imerode(logical(Ps),diamond_kernel);
Ps_edge=Ps.*(1-Ps_erode);
[edge_row,edge_col]=pupilEdge_cal(pupil,Ps,Ps_edge);
space_obj=zeros(nObj);

% Stepsize scheduling
for iIter=1:Opts.maxIter
    if contains(Opts.step_mode,'monotone')
        step_pupil(iIter)=step_max_pupil.*(Opts.maxIter-iIter+1)/(Opts.maxIter);
        step_obj(iIter)=step_max_obj.*(Opts.maxIter-iIter+1)/(Opts.maxIter);
    elseif contains(Opts.step_mode,'constant')
        step_pupil(iIter)=step_max_pupil;
        step_obj(iIter)=step_max_obj;
    end

end

loss_cal_mode=Opts.loss_cal_mode;

%% Reconstruction
while iter < Opts.maxIter
    iter = iter+1;
    loss_global_obj=(zeros(nObj));
    loss_global_pupil=(zeros(nSampR));
    for iImg=1:nImg

        %% Forward model
        I_mea = I(:,:,iImg); %measured intensity
        kyl=cropR(iImg,1);kyh=cropR(iImg,2); %corresponding x,y coordinate for each measurement
        kxl=cropR(iImg,3);kxh=cropR(iImg,4);
        highResFT_crop=highResFT(kyl:kyh,kxl:kxh);
        lowResFT=highResFT_crop.*pupil;
        lowRes = IFT(lowResFT);
        lowResFT_update = FT(sqrt(I_mea).*lowRes./(abs(lowRes+eps)));

        %% Amplitude Cost Function
        if contains(Opts.loss_mode,'amplitude')
            % obj loss
            if loss_cal_mode=='GS'
                loss_local_obj=conj(pupil).*(lowResFT_update - lowResFT) ./ (max(max((abs(pupil)).^2)));
                loss_local_pupil=conj(highResFT_crop).*(lowResFT_update - lowResFT) ./ (max(max((abs(highResFT_crop)).^2))).* Ps;
            end
                
            if loss_cal_mode=='GN'
                Opts.reg_alpha=5;
                Opts.reg_beta=5;
                loss_local_obj=abs(pupil) .* conj(pupil) .* (lowResFT_update - lowResFT) ./ max(max(abs(pupil))) ./ (abs(pupil).^2 + Opts.reg_alpha);
                loss_local_pupil=abs(highResFT_crop) .* conj(highResFT_crop) .* (lowResFT_update - lowResFT) / max(max(abs(highResFT_crop))) ./ (abs(highResFT_crop).^2 + Opts.reg_beta) .* Ps;
            end
            
            loss_global_obj(kyl:kyh,kxl:kxh)=loss_global_obj(kyl:kyh,kxl:kxh)+loss_local_obj;
            loss_global_pupil=loss_global_pupil+loss_local_pupil;
        
        elseif contains(Opts.loss_mode,'intensity')
            if Opts.adaptiveStep==1
                factor=sum(sum(I(:,:,1)))./sum(sum(I_mea(:)));
            else
                factor = 1;
            end
            %% Intensity Cost Function
            % obj loss
            loss_local_obj=FT((conj(lowRes).*lowRes-I_mea).*lowRes).*conj(pupil);
            loss_local_obj=loss_local_obj.*(factor);
            loss_global_obj(kyl:kyh,kxl:kxh)=loss_global_obj(kyl:kyh,kxl:kxh)+loss_local_obj./nImg;
    
            % pupil loss
            loss_local_pupil=FT((conj(lowRes).*lowRes-I_mea).*lowRes).*conj(highResFT_crop);
            loss_local_pupil=loss_local_pupil.*(factor);
            loss_global_pupil=loss_global_pupil+loss_local_pupil./nImg;
        end

        % Rejecting vignetting image
        % if((illumNA(iImg)<NA*Opts.vignettingFactor)&&(illumNA(iImg)>Opts.NA_obj)) % updates for vignetting measurement (reject)
        %     loss_local_obj=0;
        %     loss_local_pupil=0;
        % end

        %% local update
        if (contains(Opts.update_pupil,'local'))
            pupil = pupil+step_pupil(iter).*(loss_local_pupil).*Ps;
        end

        if (contains(Opts.update_obj,'local'))
            highResFT(kyl:kyh,kxl:kxh) = highResFT(kyl:kyh,kxl:kxh)+step_obj(iter).*(loss_local_obj);
        end


        %% solution space for data fidelity
        if iter==1
            space_obj(kyl:kyh,kxl:kxh)=space_obj(kyl:kyh,kxl:kxh)|Ps;
        end


    end


    %% global updates
    if (contains(Opts.update_obj,'global'))
        highResFT=highResFT-step_obj(iter).*((loss_global_obj));
    end
    if Opts.weight_regular_obj>0
        % regularization
        for iSub=1:Opts.nSubIter_obj
            divV = get_divergence(V0, V1);
            divW = get_divergence(W0, W1);
            highResFT=highResFT-step_obj(iter).*(weight_split_obj.*FT(-(get_sum_neighbour6(U)-4.*U)+(divV - divW)));
            U=IFT(highResFT);
            [U_0,U_1]=update_gradient_U(U);
            [V0,V1]=update_V(U_0,U_1,W0,W1,weight_regular_obj,weight_split_obj);
            [W0,W1]=update_W(U_0,U_1,W0,W1,V0,V1);
        end
    end
    
    if (contains(Opts.update_pupil,'global'))
        pupil=pupil-step_pupil(iter).*(loss_global_pupil).*Ps;
    end
    if Opts.reg_pupil==1
        for iSub=1:Opts.nSubIter_pupil
            laplacian_pupil=get_sum_neighbour6_pupil(pupil,Ps_erode,edge_row,edge_col);
            pupil=pupil-step_pupil(iter).*(-weight_regular_pupil.*laplacian_pupil).*Ps;
        end
    end
    %% calculate energy
    
    if Opts.display == 1
        display_results_simulation(highResFT,pupil,freqUV,I(:,:,1),Opts);
    end
    
        
    end
    highRes=IFT(highResFT);
end


function [U_0,U_1]=update_gradient_U(U)
U_0 = get_derivative_forward(U, 1);
U_1 = get_derivative_forward(U, 2);
end

function divK=get_divergence(K_0,K_1)
divK    = 0;
divK    = divK + get_derivative_backward(K_0, 1);
divK    = divK + get_derivative_backward(K_1, 2);
%         divK    = divK + get_derivative_backward(K_2, 3);
end

function DU=get_derivative_forward(U, axis)
Up  = shift(U, axis,'forward');
DU  = Up - U;
end



function DU=get_derivative_backward(U, axis)
Un  = shift(U, axis,'backward');
DU  = U - Un;
U_size=size(U);

if axis == 1
    DU(U_size(1), :, :) = - Un(U_size(1), :, :);

elseif axis == 2
    DU(:, U_size(2), :) = - Un(:, U_size(2), :);

elseif axis == 3
    DU(:, :, U_size(3)) = - Un(:,: , U_size(3));

end
end

function U_shift=shift(U, axis,opt)
U_size=size(U);
if contains(opt,'backward')
    if axis == 1
        shift_arr = [1, 0, 0];
    elseif axis == 2
        shift_arr = [0, 1, 0];
    elseif axis == 3
        shift_arr = [0, 0, 1];
    end
    U_shift=circshift(U,shift_arr);

    if axis == 1
        U_shift(1, :, :) = 0;
    elseif axis == 2
        U_shift(:, 1, :) = 0;
    elseif axis == 3
        U_shift(:, :, 1) = 0;

    end
elseif contains(opt,'forward')
    if axis == 1
        shift_arr = [-1, 0, 0];
    elseif axis == 2
        shift_arr = [0, -1, 0];
    elseif axis == 3
        shift_arr = [0, 0, -1];
    end
    U_shift=circshift(U,shift_arr);

    if axis == 1
        U_shift(U_size(1), :, :) = U(U_size(1), :, :);
    elseif axis == 2
        U_shift(:, U_size(2), :) = U(:, U_size(2), :);
    elseif axis == 3
        U_shift(:, :, U_size(3)) = U(:, :, U_size(3));
    end


end


end


function [pupilEdge_row,pupilEdge_col]=pupilEdge_cal(P,mask,mask_edge)
P_size=size(P);
[rowRange,colRange]=find(mask_edge);
pupilEdge_row=zeros(P_size(1),P_size(2),4);
pupilEdge_col=zeros(P_size(1),P_size(2),4);

for i = 1:length(rowRange)
    %% upleft
    if((mask(rowRange(i)-1,colRange(i))==0)&&((mask(rowRange(i)+1,colRange(i))==0)))
        pupilEdge_row(rowRange(i),colRange(i),1)=1;
    elseif((mask(rowRange(i)-1,colRange(i))==0)&&(mask(rowRange(i)+1,colRange(i))~=0))
        pupilEdge_row(rowRange(i),colRange(i),2)=1;
    elseif((mask(rowRange(i)+1,colRange(i))==0)&&(mask(rowRange(i)-1,colRange(i))~=0))
        pupilEdge_row(rowRange(i),colRange(i),3)=1;
    else
        pupilEdge_row(rowRange(i),colRange(i),4)=1;
    end
    if ((mask(rowRange(i),colRange(i)+1)==0)&&(mask(rowRange(i),colRange(i)-1)==0))
        pupilEdge_col(rowRange(i),colRange(i),1)=1;
    elseif ((mask(rowRange(i),colRange(i)-1)==0))&&(mask(rowRange(i),colRange(i)+1)~=0)
        pupilEdge_col(rowRange(i),colRange(i),2)=1;
    elseif((mask(rowRange(i),colRange(i)+1)==0)&& (mask(rowRange(i),colRange(i)-1)~=0))
        pupilEdge_col(rowRange(i),colRange(i),3)=1;
    else
        pupilEdge_col(rowRange(i),colRange(i),4)=1;
    end
end

end


function Psum=get_sum_neighbour6_pupil(P,mask_erode,pupilEdge_row,pupilEdge_col)

P_size=size(P);

%% up
Upoo    = circshift(P, [-1,0,0]);
Upoo(P_size(1),:,:)=P(P_size(1),:,:);

%% down
Unoo    = circshift(P, [1,0,0]);
Unoo(1,:,:)=P(1,:,:);

%% left
Uopo    = circshift(P, [0,-1,0]);
Uopo(:,P_size(2),:)=P(:,P_size(2),:);

%% right
Uono    = circshift(P, [0,1,0]);
Uono(:,1,:)=P(:,1,:);

Psum_row=(Upoo-P).*pupilEdge_row(:,:,2) +(Unoo-P).*pupilEdge_row(:,:,3) +(Upoo+Unoo-2.*P).*pupilEdge_row(:,:,4);
Psum_col=(Uopo-P).*pupilEdge_col(:,:,2) +(Uono-P).*pupilEdge_col(:,:,3) +(Uopo+Uono-2.*P).*pupilEdge_col(:,:,4);


Psum = (Upoo + Unoo + Uopo + Uono-4.*P).*mask_erode + Psum_row + Psum_col;


end

function Usum=get_sum_neighbour6(U)

U_size=size(U);

Upoo    = circshift(U, [-1,0,0]);
Upoo(U_size(1),:,:)=U(U_size(1),:,:);

Unoo    = circshift(U, [1,0,0]);
Unoo(1,:,:)=U(1,:,:);

Uopo    = circshift(U, [0,-1,0]);
Uopo(:,U_size(2),:)=U(:,U_size(2),:);

Uono    = circshift(U, [0,1,0]);
Uono(:,1,:)=U(:,1,:);

Usum = Upoo + Unoo + Uopo + Uono;

end


function [V0,V1,V2]=update_V(U_0,U_1,W0,W1,weight_regular,weight_split)

V0 = wthresh(real(U_0 + W0),'s',weight_regular / weight_split)+1i.*wthresh(imag(U_0 + W0),'s',weight_regular / weight_split);
V1 = wthresh(real(U_1 + W1),'s',weight_regular / weight_split)+1i.*wthresh(imag(U_1 + W1),'s',weight_regular / weight_split);
% V0 = wthresh(abs(U_0 + W0),'s',weight_regular_abs / weight_split).*exp(1i.*wthresh(angle(U_0 + W0),'s',weight_regular_angle / weight_split));
% V1 = wthresh(abs(U_1 + W1),'s',weight_regular_abs / weight_split).*exp(1i.*wthresh(angle(U_1 + W1),'s',weight_regular_angle / weight_split));



%         V2 = wthresh(U_2 + W2,'s',weight_regular / weight_split);

end

function [W0,W1,W2]=update_W(U_0,U_1,W0,W1,V0,V1)

W0 = W0 + U_0 - V0;
W1 = W1 + U_1 - V1;
%         W2 = W2 + U_2 - V2;
end


function [energy]=compute_energy(I,U,U_0,U_1,V0,V1,W0,W1,weight_regular,weight_split,FT,IFT,cropR,pupil)
dim=size(I);
highResFT=FT(U);
ampEstimate=I;
for iImg = 1:dim(3)
    kyl=cropR(iImg,1);kyh=cropR(iImg,2); %corresponding x,y coordinate for each measurement
    kxl=cropR(iImg,3);kxh=cropR(iImg,4);
    highResFT_crop=highResFT(kyl:kyh,kxl:kxh);
    lowResFT=highResFT_crop.*pupil;
    lowRes = IFT(lowResFT);
    ampEstimate(:,:,iImg)=(abs(lowRes));
end
diff=(ampEstimate - sqrt(I));
energy_data      = 0.5 * get_norm_L22(diff)/ numel(I);
energy_regular   = weight_regular * get_norm_L11([V0, V1]) / numel(I);
energy_split     = 0.5 * weight_split * get_norm_L22([U_0-V0+W0, U_1-V1+W1]) / numel(I);
energy=[energy_data + energy_regular + energy_split,energy_data,energy_regular,energy_split] ;

end

function norm_total= get_norm_L22(arr)
arr_flatten=reshape(arr,1,[]);
norm_total=sqrt(sum(abs(arr_flatten).^2));

end

function norm_total= get_norm_L11(arr)
arr_flatten=reshape(arr,1,[]);
norm_total=(sum(abs(arr_flatten)));
end