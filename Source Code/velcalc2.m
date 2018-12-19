function [vx,vy,v,X,Y,T,good,G,S,vgood]= velcalc2(XYT,resthresh,mspersamp,MINV,MAXV,xwin,ywin,twin,how_many2)

% function [vx,vy,v,x,y,t,good,G,S]=velocitycalc(XYT,resthresh,Outliers,mspersamp,MAXV,xwin,ywin,twin,how_many2)
% Use local fits to get velocity, drop points where estimates are bad
% Parameters
% XYT: 		data set from localfit which includes points and fits
% resthresh: 	threshold for residual error
% outliers:	switch to set if you wish to exclude outliers
% mspersamp:	scale velocity by sampling interval in time
% MAXV:		ignore velocities larger than MAXV
% Output:
% vx,vy,v:	velocity components and magnitude
% x,y,t:	positions in x,y,t space where good estimates were found
vx=NaN*ones(length(XYT),1);
vy=vx;
X=vx;
Y=vx;
T=vx;
G=vx;
S=vx;
col=zeros(length(XYT),1);
k1=1.0;k2=3.0;

% Extract positions from XYT set and compute velocities
L = length(XYT(:,1));
for i=1:L
   X(i)=XYT(i,1);
   Y(i)=XYT(i,2);
   T(i)=XYT(i,3);
   dx=abs(XYT(:,1)-X(i));
   dy=abs(XYT(:,2)-Y(i));
   dt=abs(XYT(:,3)-T(i));
   near=find((dx<=xwin)&(dy<=ywin)&(dt<=twin));
   xytn = -XYT(near,1:3)+ones(length(near),1)*XYT(i,1:3);
   C = XYT(near,4:9);

   % problems (don't include in avg)
   bad1=(XYT(near,10)>resthresh);	% correlation too low
   bad2=(XYT(near,13)>0.95);			% linear correlation too low
   bad3=(XYT(near,12)>1000);			% ill-conditioned fit
   bad4=(XYT(near,11)<how_many2);
   good = find(~(bad1|bad2|bad3|bad4));
   G(i) = length(good);

   if ~isempty(good)
      t=xytn(good,3);
      x=xytn(good,1);
      y=xytn(good,2);
      tx = C(good,2) + 2*C(good,4).*x + C(good,6).*y;
      ty = C(good,3) + 2*C(good,5).*y + C(good,6).*x;
      linres=XYT(near(good),13);
      res=XYT(near(good),10);
      W = k1*(1./(((linres.^2).*sqrt(x.^2 + y.^2 + t.^2)+ ...
		0.1*sqrt(xwin^2 + ywin^2 + twin^2)))) + k2*(1./(res));
      Tx = sum(tx.*W)./sum(W);
      Ty = sum(ty.*W)./sum(W);
      
      Tm=(Tx^2 + Ty^2);
      if Tm>0
         % Check consistency of gradient magnitude
         t=(tx.^2 + ty.^2);
         S(i) = std(t)/mean(t);

         vx(i) = Tx/Tm;
         vy(i) = Ty/Tm;
      end
   else
      vx(i)=NaN;
      vy(i)=NaN;
   end
       if ~isempty(good)
        if mean(linres)<.85 && mean(res)<.25
         col(i)=1;
        end
       end
end
vx=vx/mspersamp;		% correct for sampling interval
vy=vy/mspersamp;
v=sqrt(vx.^2+vy.^2);		% velocity magnitude
   
   
%Problems
bad1=(isnan(vx)|isnan(vy));		% NaN errors
bad2=(v>MAXV);
bad3=(v<MINV);                            % non-physiological velocity
bad4=S>1;                      % fits disagree
bad5=G<2;                        % too few fits
bad=(bad1|bad2|bad3|bad4|bad5);

% Do this in GUI now
% if Outliers,
%    bad6=(abs(vx)>median(abs(vx))+Outliers*std(abs(vx)));
%    bad7=(abs(vy)>median(abs(vy))+Outliers*std(abs(vy)));
%    bad8=(abs(v)>median(abs(v))+Outliers*std(abs(v)));
%    
%    bad=(bad|bad6|bad7|bad8);
% end;

good = find(~bad);

vgood=find(col);
a=zeros(length(vgood),1);
for i=1:length(vgood)
  b=find(good==vgood(i));
  if ~isempty(b)
   a(i)=b;
  end
end

aa=find(a);
vgood=a(aa);
