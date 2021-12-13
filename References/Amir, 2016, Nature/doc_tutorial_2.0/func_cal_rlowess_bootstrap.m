function [y_Bootstrap,xs]=func_cal_rlowess_bootstrap(Overlap,RootJSD,NumBootstraps,span)

NumSamples=length(Overlap);

if nargin<4
    span=0.2;
end
xs=Overlap(:);
xs(isnan(xs))=[];
% [xs,ind]=sort(xs);
xs=linspace(min(xs),max(xs),500);
if NumSamples<50
   xs=linspace(min(xs),max(xs),20);
end 

y_Bootstrap=nan(NumBootstraps, length(xs));

parfor i=1:NumBootstraps
    % resample NumSamples samples
    k=randi(NumSamples,NumSamples,1);
    ChosenOverlap=Overlap(k,k); %#ok<*PFBNS>
    ChosenRJSD=RootJSD(k,k);
    OverlapVector=ChosenOverlap(:);
    RJSDVector=ChosenRJSD(:);
    
    OverlapVector(isnan(OverlapVector))=[];
    RJSDVector(isnan(RJSDVector))=[];
    
%     ys=mylowess([OverlapVector  RJSDVector],xs,span);
    ys=myrlowess([OverlapVector  RJSDVector],xs,span);

    y_Bootstrap(i,:)=ys(:)';

end % for i

end % func_cal_rlowess_bootstrap

%%%%%%%%%%%%%%%%%%%%%%
function ys=myrlowess(xy,xs,span)
    %MYLOWESS Lowess smoothing, preserving x values
    %   YS=MYLOWESS(XY,XS) returns the smoothed version of the x/y data in the
    %   two-column matrix XY, but evaluates the smooth at XS and returns the
    %   smoothed values in YS.  Any values outside the range of XY are taken to
    %   be equal to the closest values.

    if nargin<3 || isempty(span)
        span = .3;
    end

    % Sort and get smoothed version of xy data
    xy = sortrows(xy);
    x1 = xy(:,1);
    y1 = xy(:,2);
    ys1 = smooth(x1,y1,span,'rlowess');

    % Remove repeats so we can interpolate
    t = diff(x1)==0;
    x1(t)=[]; ys1(t) = [];

    % Interpolate to evaluate this at the xs values
    ys = interp1(x1,ys1,xs,'linear',NaN);

    % Some of the original points may have x values outside the range of the
    % resampled data.  Those are now NaN because we could not interpolate them.
    % Replace NaN by the closest smoothed value.  This amounts to extending the
    % smooth curve using a horizontal line.
    if any(isnan(ys))
        ys(xs<x1(1)) = ys1(1);
        ys(xs>x1(end)) = ys1(end);
    end
end % function myrlowess