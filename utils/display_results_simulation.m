function display_results_simulation(objFT,pupil,freqUV,lowResRaw,Opts)

% Display results
if Opts.display
    objR = Opts.IFT(objFT);
    if strcmp(Opts.mode,'real')
        o = objR;
    elseif strcmp(Opts.mode,'fourier')
        o = Opts.FT(objR);
    end
    
    maximg_amp=max(max(abs(Opts.objGT)));
    minimg_amp=min(min(abs(Opts.objGT)));
    maximg_phase=max(max(angle(Opts.objGT)));
    minimg_phase=min(min(angle(Opts.objGT)));
    max_pupil=max(max(angle(Opts.pupilGT)));
    min_pupil=min(min(angle(Opts.pupilGT)));
    minimg_amp=0;
    figure(88);set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    subplot(2,5,1); imagesc(abs(Opts.objGT)); colormap gray; axis image off;  title('GT amp','FontSize',20);caxis([minimg_amp,maximg_amp]);
    subplot(2,5,2); imagesc(abs(objR)); colormap gray; axis image off;title('recon amp','FontSize',20);caxis([minimg_amp,maximg_amp]);
    subplot(2,5,3); imagesc(lowResRaw); axis image off;colormap gray; title('LR img','FontSize',20);
    

    subplot(2,5,4); imagesc(abs(pupil)); axis image off;colormap gray; title('pupil abs','FontSize',20);
    subplot(2,5,5); imagesc(log(abs(objFT))); axis image off;colormap gray; title('highRes FT','FontSize',20);

    subplot(2,5,8); imagesc(angle(Opts.pupilGT)); axis image off;colormap gray; title('pupil GT phase','FontSize',20); caxis ([min_pupil,max_pupil]);
    subplot(2,5,6); imagesc(angle(Opts.objGT)); colormap gray; axis image off;caxis([minimg_phase,maximg_phase]);title('GT phase','FontSize',20);
    subplot(2,5,7); imagesc(angle(objR)); colormap gray; axis image off;title('recon phase','FontSize',20); caxis([minimg_phase,maximg_phase]);
    subplot(2,5,9); imagesc(angle(pupil)); axis image off;colormap gray; title('pupil phase','FontSize',20); caxis([min_pupil,max_pupil]);
    
    drawnow;
end
end