function [ifr, t, rate_tl] = plotInstFR(voltage, time, transient, tl, Xlabel, colors, FigPosition)

    if nargin < 4;
        tl = xlimits;
    end

    if nargin < 5;
        Xlabel = [];
    else
        Xlabel = Xlabel;
    end

    if nargin < 6;
    colors = [
            %254/255    153/255   82/255     % naranja
            %47/255     191/255   246/255    % azul-clarito
            33/255    113/255   181/255    % azul equal instFig
            239/255   59/255    44/255     % rojo equal instFig
            ];
    else
        colors = colors;
    end


    if nargin < 7;
        FigPosition = [680 678 560 420];
    else
        FigPosition = FigPosition;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    npopulations = numel(voltage);
    for ipop = 1:npopulations
        [iV, iNeuron] = size(voltage{npopulations});
        raster{ipop} = computeRaster(time,voltage{ipop});
        pool{ipop} = 1:iNeuron;
    end

    xlimits = [time(1)+transient/2, time(end)]; % [tOn, tOff];
    uscaling      = 1e3; % unit scaling from kernel regression from kHz to Hz
    kwidth        = 5; % width of kernel regression in ms
    dt            = time(2)-time(1);
    Ts            = 1; % subsampling period in ms
    flag_interp   = 0;
    kernel        = 'L'; % 'E'; % 'G'; % for Laplacian, Epanechnikov (default) and Gaussian
    ifr           = [];


    color_line     = [178 178 178]./255; % Gris
    lineWidth     = 1.5;
    lineWidth_gfc = .2;
    fontSize      = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if npopulations == 2;
    if ~all(cellfun(@isempty,raster))
        raster = raster(~cellfun(@isempty,raster));
        rate_tl = plotRaster(raster,pool,xlimits,tl,Xlabel,colors);
        for ipop = 1:npopulations
            [ifr(:,ipop), t(:,ipop)] = NWKraster(time, raster{ipop}, pool{ipop}, kwidth, Ts, flag_interp, kernel);
        end
        ifr = uscaling*ifr;
    else
        disp('empty rasters')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    hold on
    set(gca,'layer','top')
    yline = linspace(0,100,8);
    xline = ones(1,8)*(tl(1));

    for ipop = 1:npopulations
        set(0,'DefaultFigurePosition', FigPosition)
        plot(t(:,ipop),ifr(:,ipop),'color',colors(ipop,:),'lineWidth',lineWidth)
        hold on
        plot(xline,yline,'--','color',color_line,'lineWidth', lineWidth_gfc)
    end
    xlim([transient/2 max(time)])
    ylim([0 100])
    set(gca,'fontSize',fontSize,'LineWidth',lineWidth_gfc,'TickDir','out','Box','on','XTickLabel',Xlabel,'YTick',[0:50:100]);
    xlabel('Time (s)','fontSize',fontSize)
    ylabel('iFR (sp/s)','fontSize',fontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else npopulations == 4;
        if ~all(cellfun(@isempty,raster))
            raster = raster(~cellfun(@isempty,raster));
            rate_tl{1} = plotRaster(raster(1:2),pool(1:2),xlimits,tl,Xlabel,colors);
            rate_tl{2} = plotRaster(raster(3:4),pool(3:4),xlimits,tl,Xlabel,colors);
            for ipop = 1:npopulations
                [ifr(:,ipop), t(:,ipop)] = NWKraster(time, raster{ipop}, pool{ipop}, kwidth, Ts, flag_interp, kernel);
            end
            ifr = uscaling*ifr;
        else
            disp('empty rasters')
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    hold on
    set(gca,'layer','top')
    yline = linspace(0,100,8);
    xline = ones(1,8)*(tl(1));
                for ipop = 1:2
            set(0,'DefaultFigurePosition', FigPosition)
            plot(t(:,ipop),ifr(:,ipop),'color',colors(ipop,:),'lineWidth',lineWidth)
            hold on
            plot(xline,yline,'--','color',color_line,'lineWidth', lineWidth_gfc)
        end
    xlim([transient/2 max(time)])
    ylim([0 100])
    set(gca,'fontSize',fontSize,'LineWidth',lineWidth_gfc,'TickDir','out','Box','on','XTickLabel',Xlabel,'YTick',[0:50:100]);
    xlabel('Time (s)','fontSize',fontSize)
    ylabel('iFR (sp/s)','fontSize',fontSize)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    hold on
    set(gca,'layer','top')
        for ipop = 3:4
            set(0,'DefaultFigurePosition', FigPosition)
            plot(t(:,ipop),ifr(:,ipop),'color',colors(ipop-2,:),'lineWidth',lineWidth)
            hold on
            plot(xline,yline,'--','color',color_line,'lineWidth', lineWidth_gfc)
        end
    xlim([transient/2 max(time)])
    ylim([0 100])
    set(gca,'fontSize',fontSize,'LineWidth',lineWidth_gfc,'TickDir','out','Box','on','XTickLabel',Xlabel,'YTick',[0:50:100]);
    xlabel('Time (s)','fontSize',fontSize)
    ylabel('iFR (sp/s)','fontSize',fontSize)
    end

