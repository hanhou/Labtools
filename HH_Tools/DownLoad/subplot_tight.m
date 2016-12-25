function h=subplot_tight(m,n,p,gaps,margins,varargin)
%function subplot_tight(m,n,p,gaps,margins,varargin)
%
% Functional purpose: A wrapper function for Matlab function subplot. Adds the ability to define the gaps between
% neighbouring subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the gaps between
% subplots can reach 40% of figure area, which is pretty lavish.  
%
% Input arguments (defaults exist):
%   gaps- two elements vector [vertical,horizontal] defining the gaps between neighbouring axes. Default value
%            is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%            relatively large axis. 
%   margins - [UP DOWN LEFT RIGHT]
%
% Output arguments: same as subplot- none, or axes handle according to function call.
%
% Issues & Comments: Note that if additional elements are used in order to be passed to subplot, gaps parameter must
%       be defined. For default gaps value use empty element- [].      
%
% Author and Date:  Nikolay S. 29/03/2011. 
% Last update:      Nikolay S. 21/04/2011 (accourding to Alan B comment).
%
% Usage example: h=subplot_tight((2,3,1:2,[0.5,0.2])

if (nargin<4) || isempty(gaps)
    gaps=[0.01,0.01]; % default gaps value- 1% of figure
end

if length(gaps)==1
    gaps(2)=gaps;
end

if (nargin<5) || isempty(margins) 
    margins = [gaps(1) gaps(1) gaps(2) gaps(2)];
elseif length(margins) < 4
    margins = repmat(margins(1),1,4);
end


%note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
[subplot_col,subplot_row]=ind2sub([n,m],p);  


height=(1-(m-1)*gaps(1)-margins(1)-margins(2))/m; % single subplot height
width=(1-(n-1)*gaps(2)-margins(3)-margins(4))/n;  % single subplot width

% note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot 
subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot   

merged_height=subplot_rows*( height+gaps(1) )- gaps(1);   % merged subplot height
merged_width= subplot_cols*( width +gaps(2) )- gaps(2);   % merged subplot width

merged_bottom=(m-max(subplot_row))*(height+gaps(1)) + margins(2); % merged subplot bottom position
merged_left=(min(subplot_col)-1)*(width+gaps(2)) + margins(3);              % merged subplot left position
pos_vec=[merged_left merged_bottom merged_width merged_height];

% h_subplot=subplot(m,n,p,varargin{:},'Position',pos_vec);
% Above line doesn't work as subplot tends to ignore 'position' when same mnp is utilized
h_subplot=subplot('Position',pos_vec,varargin{:});

if nargout~=0
    h=h_subplot;
end