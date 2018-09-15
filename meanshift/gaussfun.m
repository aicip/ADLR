function out = gaussfun(x,d,bandWidth)
% approximate Gaussian kernel
% x - value
% d - location
% bandWidth - band width of the kernel
% out - contribtion to the kernel mean
%
% Copyright 2015 Han Gong, University of East Anglia
%     persistent ns = 1000; % resolution of guassian approximation
%     persistent xs = linspace(0,bandWidth,ns+1); % approximate ticks
%     persistent kfun = exp(-(xs.^2)/(2*bandWidth^2));
    ns = 1000; % resolution of guassian approximation
    xs = linspace(0,bandWidth,ns+1); % approximate ticks
    kfun = exp(-(xs.^2)/(2*bandWidth^2));
    w = kfun(round(d/bandWidth*ns)+1);
    w = w/sum(w); % normalize
    out = sum( bsxfun(@times, x, w ), 2 );
end
