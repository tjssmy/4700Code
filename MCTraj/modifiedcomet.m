function modifiedcomet(varargin)
%MODIFIEDCOMET  Comet-like trajectory, modified from original comet.m

%   COMET(Y) displays an animated comet plot of the vector Y.
%   COMET(X,Y) displays an animated comet plot of vector Y vs. X.
%   COMET(X,Y,p) uses a comet of length p*length(Y).  Default is p = 0.10.

%   COMET(X,Y,p,num) uses a comet of length p*length(Y), num determines the
%   number of trajectories to plot. Max is 5.

%   COMET(AX,...) plots into AX instead of GCA.
%
%   Example:
%       t = -pi:pi/200:pi;
%       comet(t,tan(sin(t))-sin(tan(t)))
%
%   See also COMET3.

%   Charles R. Denham, MathWorks, 1989.
%   Revised 2-9-92, LS and DTP; 8-18-92, 11-30-92 CBM.
%   Copyright 1984-2015 MathWorks, Inc.

% Parse possible Axes input
[ax,args,nargs] = axescheck(varargin{:});
if nargs < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nargs > 4 % modified to allow for 4 inputs. New input is "num".
    error(message('MATLAB:narginchk:tooManyInputs'));
end

% Parse the rest of the inputs
if nargs < 2, x = args{1}; y = x; x = 1:length(y); end
if nargs == 2, [x,y] = deal(args{:}); end
if nargs < 3, p = 0.10; end
if nargs == 3, [x,y,p] = deal(args{:}); end
if nargs == 4, [x,y,p,num] = deal(args{:}); end % modified to include "num"

if ~isscalar(p) || ~isreal(p) ||  p < 0 || p >= 1
    error(message('MATLAB:comet:InvalidP'));
end

ax = newplot(ax);
if ~ishold(ax)
    [minx,maxx] = minmax(x);
    [miny,maxy] = minmax(y);
    axis(ax,[minx maxx miny maxy])
end
co = get(ax,'colororder');

m = length(x);
k = round(p*m);

if size(co,1)>=3
    if num == 1
        colors = [ co(1,:);co(6,:);co(1,:)]; %First trajectory color
    end
    if num == 2
        colors = [ co(2,:);co(6,:);co(2,:)]; %Second trajectory color
    end
    if num == 3
        colors = [ co(3,:);co(6,:);co(3,:)]; %Third trajectory color
    end
    if num == 4
        colors = [ co(4,:);co(6,:);co(4,:)]; %Fourth trajectory color
    end
    if num == 5
        colors = [ co(5,:);co(6,:);co(5,:)]; %Fifth trajectory color
    end
    
    lstyle = '-';
else
    colors = repmat(co(1,:),3,1);
    lstyle = '--';
end

head = line('parent',ax,'color',colors(1,:),'marker','o','linestyle','none', ...
            'xdata',x(1),'ydata',y(1),'Tag','head');

body = matlab.graphics.animation.AnimatedLine('color',colors(2,:),...
    'linestyle',lstyle,...
    'Parent',ax,...
    'MaximumNumPoints',max(1,k),'tag','body');
tail = matlab.graphics.animation.AnimatedLine('color',colors(3,:),...
    'linestyle','-',...
    'Parent',ax,...
    'MaximumNumPoints',1+m,'tag','tail'); %Add 1 for any extra points

if ( length(x) < 2000 )
    updateFcn = @()drawnow;
else
    updateFcn = @()drawnow('update');
end

% Grow the body
for i = 1:k
    set(head,'xdata',x(i),'ydata',y(i));
    if  ~( body.isvalid() )
        return
    end
    addpoints(body,x(i),y(i));
    updateFcn();

end
% Add a drawnow to capture any events / callbacks
drawnow;
% Primary loop
for i = k+1:m
    set(head,'xdata',x(i),'ydata',y(i));
    if ~( body.isvalid() )
        return
    end
    addpoints(body,x(i),y(i));
    addpoints(tail,x(i-k),y(i-k));
    updateFcn();
end
drawnow;
% Clean up the tail
for i = m+1:m+k
    if  ~( body.isvalid() )
        return
    end
    addpoints(tail,x(i-k),y(i-k));
    updateFcn();
end
drawnow;
end

function [minx,maxx] = minmax(x)
minx = min(x(isfinite(x)));
maxx = max(x(isfinite(x)));
if minx == maxx
    minx = maxx-1;
    maxx = maxx+1;
end
end