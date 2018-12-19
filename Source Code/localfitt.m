function XYT = localfitt(x,y,t,winx,winy,wint,how_many,MINt,MAXt)

% function XYT=localfit(x,y,t,winx,winy,wint,how_many,MINt,MAXt)
% Use local fits to get velocity
% Parameters
% x,y,t:	locations of activity in space and time
% win: 	size of region to fit (x, y, and t)
% how_many: 	number of points needed for good fit (use # of coefficients of unknown
% MINt, MAXt: 	start and end times
% Output:
% XYT:		matrix containing coefficients of fits and statistics of fit

alpha = 0.0;

xyt=[x,y,t];
XYT=zeros(length(xyt),14);
chk = zeros(length(XYT),1);

% Generate list of points to be fitted
fitlist=find((xyt(:,3)>=MINt)&(xyt(:,3)<=MAXt));
M=length(fitlist);

% Main algorithm: fit
for j=1:M
    i=fitlist(j);
    dx=abs(xyt(:,1)-xyt(i,1));
    dy=abs(xyt(:,2)-xyt(i,2));
    dt=abs(xyt(:,3)-xyt(i,3));
    
    near=find((dx<=winx)&(dy<=winy)&(dt<=wint));
    len = length(near);
    if len>how_many
        xytn = xyt(near,:)-ones(len,1)*xyt(i,:);
        t=xytn(:,3);
        x=xytn(:,1);
        y=xytn(:,2);
        fit = [ones(len,1) x y (alpha+x.^2) (alpha+y.^2) (alpha+x.*y)];	%6 coefficients
        coefs = fit\t;
        var_t = sum((t-mean(t)).^2);
        resi    = sqrt(sum((t-fit*coefs).^2)/var_t);
        resilin = sqrt(sum((t-fit(:,1:3)*coefs(1:3)).^2)/var_t);
        XYT(i,:)=[xyt(i,:),coefs',resi,len,cond(fit),resilin,sqrt(var_t)];
        chk(i)=1;
    end
end

XYT=XYT(find(chk),:);


%XYT--1-3 xyt
%     4-9 coefficients
%     10 resi(residual error)
%     11 len (number of points in fit
%     12 cond(fit) condition of fit
%     13 resilin (linear residual error)
%     14 var_t variance with resprect to t







