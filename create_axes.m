function [ax] = create_axes(npanels,orientation,varargin)

params.lowercase = false;
params.color = 'k';
for i = 1:2:length(varargin) % work for a list of name-value pairs
    if ischar(varargin{i}) % check if is character
        params.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
    end
end

if params.lowercase
    labels = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o',...
          'p','q','r','s','t','u','v','w','x','y','z'];
else
    labels = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',...
          'P','Q','R','S','T','U','V','W','X','Y','Z'];
end

switch orientation
    case 'vertical'
    ax = gobjects(1,npanels);
    h = (0.85-0.05*(npanels-1))/npanels;
    for npanel = 1:npanels
        ax(npanel) = axes('Position',[0.1 0.1+(h + 0.05)*(npanels-npanel) 0.8 h]);
        annotation('textbox',[0.11,(h + 0.05)*npanel,0.05,0.05],'String',labels(npanels-npanel+1),'FaceAlpha',0,'LineStyle','none','FontSize',18,'Color',params.color)
    end

    case 'double_column'
    ncols = 2;
    nrows = ceil(npanels/ncols);
    ax = gobjects(nrows,ncols);

    w = (0.85-0.05*(ncols-1))/ncols;
    h = (0.85-0.05*(nrows-1))/nrows;
    for npanel = 1:npanels
        nrow = floor((npanel-1)/ncols) + 1;
        ncol = rem(npanel-1,ncols)+1;
            ax(nrow,ncol) = axes('Position',[0.1+(w+0.025)*(ncol-1)-0.025 0.1+(h + 0.05)*(nrows-nrow) w h]);
            annotation('textbox',[0.11+(w+0.025)*(ncol-1)-0.025,(h + 0.05)*(nrows-nrow+1),0.05,0.05],...
                'String',labels(npanel),'FaceAlpha',0,'LineStyle','none','FontSize',18,'Color',params.color)
    end
    ax = ax';

    case 'double_row'
    nrows = 2;
    ncols = ceil(npanels/nrows);
    ax = gobjects(nrows,ncols);

    w = (0.85-0.05*(ncols-1))/ncols;
    h = (0.85-0.05*(nrows-1))/nrows;
    for npanel = 1:npanels
        nrow = floor((npanel-1)/ncols) + 1;
        ncol = rem(npanel-1,ncols)+1;
            ax(nrow,ncol) = axes('Position',[0.1+(w+0.025)*(ncol-1)-0.025 0.1+(h + 0.05)*(nrows-nrow) w h]);
            annotation('textbox',[0.11+(w+0.025)*(ncol-1)-0.025,(h + 0.05)*(nrows-nrow+1),0.05,0.05],...
                'String',labels(npanel),'FaceAlpha',0,'LineStyle','none','FontSize',18,'Color','k','Color',params.color)
    end
    ax = ax';

    case 'horizontal'
    ax = gobjects(npanels,1);
    w = (0.85-0.05*(npanels-1))/npanels;
    for npanel = 1:npanels
        ax(npanel) = axes('Position',[0.1+(w + 0.05)*(npanel-1) 0.1 w 0.8]);
        annotation('textbox',[0.11+(w + 0.05)*(npanel-1) ,0.85,0.05,0.05],'String',labels(npanel),'FaceAlpha',0,'LineStyle','none','FontSize',18,'Color',params.color)
    end

    case 'horizontal_color'
    ax = gobjects(npanels,1);
    w = (0.85-0.05*(npanels-1))/npanels;
    for npanel = 1:npanels
        ax(npanel) = axes('Position',[0.05+(w + 0.05)*(npanel-1) 0.1 w 0.8]);
        annotation('textbox',[0.06+(w + 0.05)*(npanel-1) ,0.85,0.05,0.05],'String',labels(npanel),'FaceAlpha',0,'LineStyle','none','FontSize',18,'Color',params.color)
    end

    case 'tile_close'
    ncols = floor(sqrt(npanels));
    nrows = ceil(npanels/ncols);
    ax = gobjects(nrows,ncols);

    w = (0.85-0.05*(ncols-1))/ncols;
    h = (0.85-0.05*(nrows-1))/nrows;
    for npanel = 1:npanels
        nrow = floor((npanel-1)/ncols) + 1;
        ncol = rem(npanel-1,ncols)+1;
            ax(nrow,ncol) = axes('Position',[0.1+(w+0.05)*(ncol-1) 0.1+(h + 0.05)*(nrows-nrow) w h]);
            annotation('textbox',[0.11+(w+0.05)*(ncol-1),(h + 0.05)*(nrows-nrow+1),0.05,0.05],...
                'String',labels(npanel),'FaceAlpha',0,'LineStyle','none','FontSize',18,'Color',params.color)
    end
    ax = ax';

    case 'tile_close_color'
    ncols = floor(sqrt(npanels));
    nrows = ceil(npanels/ncols);
    ax = gobjects(nrows,ncols);

    w = (0.79-0.05*(ncols-1))/ncols;
    h = (0.85-0.05*(nrows-1))/nrows;
    for npanel = 1:npanels
        nrow = floor((npanel-1)/ncols) + 1;
        ncol = rem(npanel-1,ncols)+1;
            ax(nrow,ncol) = axes('Position',[0.1+(w+0.05)*(ncol-1) 0.1+(h + 0.05)*(nrows-nrow) w h]);
            annotation('textbox',[0.11+(w+0.05)*(ncol-1),(h + 0.05)*(nrows-nrow+1),0.05,0.05],...
                'String',labels(npanel),'FaceAlpha',0,'LineStyle','none','FontSize',18,'Color',params.color)
    end
    ax = ax';

    case 'tile_open'
    ncols = floor(sqrt(npanels));
    nrows = ceil(npanels/ncols);
    ax = gobjects(nrows,ncols);

    w = (0.85-0.1*(ncols-1))/ncols;
    h = (0.85-0.1*(nrows-1))/nrows;
    for npanel = 1:npanels
        nrow = floor((npanel-1)/ncols) + 1;
        ncol = rem(npanel-1,ncols)+1;
            ax(nrow,ncol) = axes('Position',[0.1+(w+0.1)*(ncol-1) 0.1+(h + 0.1)*(nrows-nrow) w h]);
            annotation('textbox',[0.11+(w+0.1)*(ncol-1),(h + 0.1)*(nrows-nrow+1)-0.05,0.05,0.05],...
                'String',labels(npanel),'FaceAlpha',0,'LineStyle','none','FontSize',18,'Color',params.color)
    end
    ax = ax';

    case 'tile_open_sharex'
    ncols = floor(sqrt(npanels));
    nrows = ceil(npanels/ncols);
    ax = gobjects(nrows,ncols);

    w = (0.85-0.1*(ncols-1))/ncols;
    h = (0.85-0.05*(nrows-1))/nrows;
    for npanel = 1:npanels
        nrow = floor((npanel-1)/ncols) + 1;
        ncol = rem(npanel-1,ncols)+1;
            ax(nrow,ncol) = axes('Position',[0.1+(w+0.1)*(ncol-1) 0.1+(h + 0.05)*(nrows-nrow) w h]);
            annotation('textbox',[0.11+(w+0.1)*(ncol-1),(h + 0.05)*(nrows-nrow+1),0.05,0.05],...
                'String',labels(npanel),'FaceAlpha',0,'LineStyle','none','FontSize',18,'Color',params.color)
    end
    ax = ax';
end

for axi = ax(1:npanels)
    set(axi,'box','on')
end
    
end