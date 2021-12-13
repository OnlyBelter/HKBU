function [xc,LocalSlope]=detect_negative_slope_Lowess_03(xx,yy)
    
% N=length(yy);
% 
% WindowsNumber=10;%100;
% WindowsLimits=linspace(min(xx),max(xx),WindowsNumber+1);


Smoothed_yy=smooth(yy);

% LocalSlope=diff(yy)./diff(xx);
LocalSlope=diff(Smoothed_yy(:))./diff(xx(:));

xc=xx(find(LocalSlope>0,1,'last'));

if isempty(xc), xc=xx(1); end


