function hSeg = segmarker(state, varargin)
%SEGMARKER SEGment MARKER

% by Gurdy Woo, 2014-12-2, 2014-12-3
% gurdy.woo@gmail.com
%
% inspired by <mhirsch@mathworks.com Michelle Hirsch>'s 
% <http://www.mathworks.com/matlabcentral/fileexchange/1542-datalabel-state-marker-color--
% datalabel>
% 

if nargin < 1
    state = 'on';
end

switch state
    case 'on'
        % Set the WindowButtonDownFcn
        set(gcf, 'WindowButtonDownFcn', @markDownCallbackFunc);
        set(gcf, 'DoubleBuffer', 'on'); % eliminate flicker    
    case 'off'
        % Unset the WindowButton{Down|Up}Fcn callback functions
        set(gcf, 'WindowButtonDownFcn', '', 'WindowButtonUpFcn', '')
    case 'make'
        if nargin == 4
            hSeg = makeSegmarker(varargin{:}); % pos, mType, pHandle
        end
    case 'delete'
        if nargin > 1
            if ischar(varargin{1})
                if strcmp(varargin{1}, 'all')
                    hSeg = findobj(gcf, '-regexp', ...
                        'Tag', 'segmarker_(foothill|valley|peak)');
                    delete(hSeg)
                end
            else
                hSeg = [varargin{:}];
                if ishandle(hSeg)
                    delete(hSeg)
                end
            end
        end
end


%-------------------------------------------------------------------------%
function tf = isSegmarker(obj)
% test a segmarker on its tag name

tagName = get(obj, 'Tag');
tf = any(strcmp(tagName, {'segmarker_foothill', 'segmarker_valley', 'segmarker_peak'}));


%-------------------------------------------------------------------------%
function tf = isLineSelected(obj)
% test a line selected excluding the segmarker

tf = strcmp(get(obj, 'Type'), 'line') && ~isSegmarker(obj);


%-------------------------------------------------------------------------%
function hSeg = makeSegmarker(pos, mType, pHandle)
% create a new segmarker
%
%  Tags:
%         tag name           mType     marker     color
%    segmarker_foothill  --    0         o          r
%    segmarker_peak      --    1         ^          r
%    segmarker_valley    --   -1         v          r
%

if nargin < 3
    pHandle = getSegline;
end

if nargin < 2
    mType = 0; % default marker type is foothill
end

tagName = 'segmarker_'; % prefix of the segmarker tag name

if mType == 0
    marker = 'o'; % circle
    tagName = [tagName 'foothill'];
elseif mType == 1
    marker = '^'; % upper triangle
    tagName = [tagName 'peak'];
elseif mType == -1
    marker = 'v'; % lower triangle
    tagName = [tagName 'valley'];
end

color = 'r'; % in red

xv = pos(1);
yv = pos(2);
hSeg = line(xv, yv, ...
    'Marker', marker, 'Color', color, 'LineStyle', 'none');

% contextmenu (right click) for basic operations
hSegmenu = uicontextmenu('Tag', 'uiSegmarkerContextmenu');

fStr = 'set(gco, ''Tag'', ''segmarker_foothill'', ''marker'', ''o'')';
uimenu(hSegmenu, 'Label', 'Change to foothill', 'Callback', fStr);
pStr = 'set(gco, ''Tag'', ''segmarker_peak'', ''marker'', ''^'')';
uimenu(hSegmenu, 'Label', 'Change to peak', 'Callback', pStr);
vStr = 'set(gco, ''Tag'', ''segmarker_valley'', ''marker'', ''v'')';
uimenu(hSegmenu, 'Label', 'Change to valley', 'Callback', vStr);

uimenu(hSegmenu, 'Label', 'Delete this segmarker', ...
    'Separator', 'on', ...
    'Callback', 'delete(gco)');
oStr = 'h1 = getappdata(gcf, ''hOutOfLineSegs''); setappdata(gcf, ''hOutOfLineSegs'', [gco h1])';
uimenu(hSegmenu, 'Label', 'Out of line', 'Callback', oStr);

set(hSeg, 'Tag', tagName, 'UserData', pHandle, 'UIContextMenu', hSegmenu);

setappdata(gcf, 'hOutOfLineSegs', [])



%-------------------------------------------------------------------------%
function markDownCallbackFunc(src, evt) %#ok
% mark down points at a line

if isLineSelected(gco)
    % User-selected point
    cp = get(gca, 'CurrentPoint');
    x = cp(1, 1);
    y = cp(1, 2);
    
    % X/Y-Data of the target line
    xl = get(gco, 'XData');
    yl = get(gco, 'YData');
    
    % Look up nearest value on the target line.
    [xv, yv] = localNearestPoint(x, xl, y, yl);
    
    makeSegmarker([xv, yv], 0, gco);

elseif isSegmarker(gco)
    
    hMenus = get(get(gco, 'UIContextMenu'), 'Children');
    
    % tricks: use strcmp
    % the uimenu handles is sorted in descending order
    % hide the one who will change to itself
    menuFlt = ~[0 0 strcmp(get(gco, 'marker'), {'v', '^', 'o'})];
    set(hMenus(menuFlt), 'Visible', 'on');
    set(hMenus(~menuFlt), 'Visible', 'off');

    % add a new menu to switch between outside and inside of the line
    hOutOfLineSegs = getappdata(gcf, 'hOutOfLineSegs');
    if ismember(gco, hOutOfLineSegs)
        if strcmp(get(hMenus(1), 'Label'), 'Out of line')
            bStr = ['h1 = getappdata(gcf, ''hOutOfLineSegs'');', ...
                '[~, b1] = ismember(gco, h1);', ...
                'h1(b1(:)) = [];', ...
                'setappdata(gcf, ''hOutOfLineSegs'', h1)'];
            
            set(hMenus(1), 'Label', 'Back to line', 'Callback', bStr);
        end
    else
        if strcmp(get(hMenus(1), 'Label'), 'Back to line')
            oStr = 'h1 = getappdata(gcf, ''hOutOfLineSegs''); setappdata(gcf, ''hOutOfLineSegs'', [gco h1])';
            set(hMenus(1), 'Label', 'Out of line', 'Callback', oStr);
        end
    end

        
    
    % unset the WindowButtonMotionFcn callback function
    btnUpCallbackFunc = @(src, evt) set(gcf, 'WindowButtonMotionFcn', '');
    
    set(gcf, 'WindowButtonMotionFcn', @dragSegmarkerCallbackFunc, ...
        'WindowButtonUpFcn', btnUpCallbackFunc);
end


%-------------------------------------------------------------------------%
function dragSegmarkerCallbackFunc(src, evt) %#ok
% drag and drop datalabels

cp = get(gca,'CurrentPoint');
x = cp(1,1);
y = cp(1,2);

% Constrain to Line
hLine = get(gco, 'UserData');

xl = get(hLine, 'XData');
yl = get(hLine, 'YData');

% Get nearest XY-value
[xv, yv] = localNearestPoint(x, xl, y, yl);

hOutOfLineSegs = getappdata(gcf, 'hOutOfLineSegs');
if ~ismember(gco, hOutOfLineSegs)
    set(gco, 'XData', xv, 'YData', yv);
else
    set(gco, 'XData', xv, 'YData', y);
end

drawnow


%-------------------------------------------------------------------------%
function [xv, yv]=localNearestPoint(x, xl, y, yl)
% Find nearest value of [xl,yl] to (x,y)
%  This function is forked from *datalabel*.
%  Inputs:
%    x   Selected x value
%    xl  Line Data (x)
%    y   Selected y value
%    yl  Line Data (y)

% Special Case: Line has a single non-singleton value
if sum(isfinite(xl))==1
    fin = find(isfinite(xl));
    xv = xl(fin);
    yv = yl(fin);
else
    % Normalize axes
    xlmin = min(xl);
    xlmax = max(xl);
    ylmin = min(yl);
    ylmax = max(yl);
    
	% Process the case where max == min
	if xlmax == xlmin
		xln = (xl - xlmin);
		xn = (x - xlmin);
	else
		% Normalize data
		xln = (xl - xlmin)./(xlmax - xlmin);
		xn = (x - xlmin)./(xlmax - xlmin);
	end
    
	if ylmax == ylmin
		yln = (yl - ylmin);
		yn = (y - ylmin);
	else
		yln = (yl - ylmin)./(ylmax - ylmin);
		yn = (y - ylmin)./(ylmax - ylmin);
	end

    % Find nearest point using our friend Ptyhagoras
    a = xln - xn;       % Distance between x and the line
    b = yln - yn;       % Distance between y and the line
    c = (a.^2 + b.^2);  % Distance between point and line
    % Don't need sqrt, since we get same answer anyway
    [~,ind] = min(c);
    
    % Nearest value on the line
    xv = xl(ind);
    yv = yl(ind);
end

