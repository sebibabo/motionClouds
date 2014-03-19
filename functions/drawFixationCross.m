function img = drawFixationCross(img, radius,value)
  % draw fixation cross
  if nargin <3 
    value = mean(img(:));
  end

  if nargin <2 
    radius = 25*sqrt(2);
  end
  
  dummy  = zeros(size(img,1),size(img,2));
  %%
  
  center = ceil([size(dummy,1) size(dummy,2)]/2);
  % white area
  % 1. draw line
  point1 = center-radius/sqrt(2);
  point2 = center+radius/sqrt(2)-1;

  rpts = linspace(point1(1),point2(1),1000);   %# A set of row points for the line
  cpts = linspace(point1(2),point2(2),1000);   %# A set of column points for the line
  index = sub2ind(size(dummy),round(rpts),round(cpts));  %# Compute a linear index
  dummy(index) = 0.5;        

  % 2. draw the other line
  point1 = center-[-1 1]*radius/sqrt(2);
  point2 = center+[-1 1]*radius/sqrt(2)-[-1 1];

  rpts = linspace(point1(1),point2(1),1000);   %# A set of row points for the line
  cpts = linspace(point1(2),point2(2),1000);   %# A set of column points for the line
  index = sub2ind(size(dummy),round(rpts),round(cpts));  %# Compute a linear index
  dummy(index) = 0.5;       

  
  % outer avg area
%   dummy  = zeros(size(img,1),size(img,2));

  point1 = center-radius/sqrt(2)+[0 -1];
  point2 = center+radius/sqrt(2)+[0 -1 ];

  rpts = linspace(point1(1),point2(1),1000);   %# A set of row points for the line
  cpts = linspace(point1(2),point2(2),1000);   %# A set of column points for the line
  index = sub2ind(size(dummy),round(rpts),round(cpts));  %# Compute a linear index
  dummy(index) = 1;        

  point1 = center-radius/sqrt(2)-[1 0];
  point2 = center+radius/sqrt(2)- [1 0];

  rpts = linspace(point1(1),point2(1),1000);   %# A set of row points for the line
  cpts = linspace(point1(2),point2(2),1000);   %# A set of column points for the line
  index = sub2ind(size(dummy),round(rpts),round(cpts));  %# Compute a linear index
  dummy(index) = 1;        
    % 2. draw the other line
  point1 = center-[-1 1]*radius/sqrt(2)+[1 0 ];
  point2 = center+[-1 1]*radius/sqrt(2)+[1 0];

  rpts = linspace(point1(1),point2(1),1000);   %# A set of row points for the line
  cpts = linspace(point1(2),point2(2),1000);   %# A set of column points for the line
  index = sub2ind(size(dummy),round(rpts),round(cpts));  %# Compute a linear index
  dummy(index) = 1;       

  point1 = center-[-1 1]*radius/sqrt(2)-[0 1];
  point2 = center+[-1 1]*radius/sqrt(2)-[0 1];

  rpts = linspace(point1(1),point2(1),1000);   %# A set of row points for the line
  cpts = linspace(point1(2),point2(2),1000);   %# A set of column points for the line
  index = sub2ind(size(dummy),round(rpts),round(cpts));  %# Compute a linear index
  dummy(index) = 1;       
  
  
  dummy = repmat(dummy,[1,1,size(img,3)]);
  img(dummy == 0.5) = 1;
  img(dummy == 1) = value;
  
  