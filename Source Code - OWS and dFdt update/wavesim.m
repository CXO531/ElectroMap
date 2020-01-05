function [p]=wavesim(wsmat,mask,ep,beforeframes);

%% Chris O'Shea, 05/09/2018
% Function for producing waveform similarty map from 4D matix of
% (row,col,frame,windowed signal)

[rows, cols, num, Nsig] = size(wsmat)
%Nsig=Nsig-1;


%% Normalise signals by standard norm
for N=1:Nsig %ignore last AP
    stnorm(:,:,:,N)=imcomplement(squeeze(wsmat(:,:,:,N)));
    %     sigsmat=imcomplement(sigsmat);
    %     sigsmat=sigsmat-mean(sigsmat,3);
    %     stnorm1=sigsmat.^2;
    %     stnorm1=sum(stnorm1,3);
    %     stnorm1=sqrt(stnorm1);
    %     stnorm(:,:,:,N)=sigsmat./stnorm1;
    
    %%shift to match
    
    %     [m,i1]=max(abs(diff(sigsmat./stnorm1,[],3)),[],3); %dv/dt
    %     %[m,i1]=min(sigsmat./stnorm1,[],3); %max
    %             for r =1:rows
    %             for c =1:cols
    %             s1=circshift(squeeze(stnorm(r,c,:,N)),-(i1(r,c)-10),1);
    %             stnorm(r,c,:,N)=s1;
    %             end
    %             end
end

% figure,
for i=1:Nsig
    %     hold on
    %     plot(squeeze(stnorm(27,17,:,i)),'b')
    %     plot(squeeze(stnorm(36,43,:,i)),'r')
    assignin('base','signals',stnorm)
end
% 
% figure,
% for i=1:Nsig
%     hold on
%     plot(squeeze(stnorm(15,15,:,i)))
% end


p=zeros(rows,cols);
for i=1:Nsig
    for j=i+1:Nsig
        sis=squeeze(stnorm(:,:,:,i));
        [min1,i1]=min(sis(:,:,1:beforeframes),[],3);
        [min12,i12]=min(sis(:,:,beforeframes:end),[],3);
        sjs=squeeze(stnorm(:,:,:,j));
        [min2,i2]=min(sjs(:,:,1:beforeframes),[],3);
        [min22,i22]=min(sjs(:,:,beforeframes:end),[],3);
        %to minimise effects of previous peak
        
        
        for r=1:rows
            for c=1:cols
                if i1(r,c)<=i2(r,c)
                    ind1(r,c)=i1(r,c);
                else
                    ind1(r,c)=i2(r,c);
                end
                
                %% code for min after as well, not used atm
                if i12(r,c)<=i22(r,c)
                    ind2(r,c)=i12(r,c)+beforeframes-1;
                else
                    ind2(r,c)=i22(r,c)+beforeframes-1;
                end
                
            end
        end
        
        
        for r=1:rows
            for c=1:cols
                %sig1=squeeze(sis(r,c,:));
                %sig2=squeeze(sjs(r,c,:));
                sig1=squeeze(sis(r,c,ind1(r,c):ind2(r,c)));
                sig2=squeeze(sjs(r,c,ind1(r,c):ind2(r,c)));
                %                 sig1=squeeze(sis(r,c,:));
                %                 sig2=squeeze(sjs(r,c,:));
                
                
                
                
                sig1=sig1-mean(sig1);
                sig11=sig1.^2;
                sig11=sum(sig11);
                sig11=sqrt(sig11);
                sig1=sig1./sig11;
                sig2=sig2-mean(sig2);
                sig21=sig2.^2;
                sig21=sum(sig21);
                sig21=sqrt(sig21);
                sig2=sig2./sig21;
                
                %% align signals by max - a bit crap beacuse of noise (increases distance even at 170)
%                 [m1,i1]=max(sig1);
%                 [m2,i2]=max(sig2);
%                 diff_max=i1-i2;
%                 sig2=circshift(sig2,diff_max);
                
                %% align signals by dv/dt
%                 [m1,i1]=max(diff(sig1));
%                 [m2,i2]=max(diff(sig2));
% %                 if isempty(i1)==1
% %                     i1
% %                     i2
% %                     ind1(r,c)
% %                     ind2(r,c)
% %                     figure,
% %                     plot(sig1b,'b')
% %                     hold on
% %                     plot(sig2b,'r')
% %                 end
%                 diff_max=i1-i2;
%                 diff_max;
%                 sig2=circshift(sig2,diff_max);
%%
                
%                 dotproduct1=(dot(sig1,sig2));
%                 
%                 d(r,c)=acos(dotproduct1);
                xy   = dot(sig1,sig2);
% %                 nx   = norm(sig1);
% %                 ny   = norm(sig2);
% nx=sig1; ny=sig2;
%                 nxny = nx*ny;
                %d(r,c)   = xy/nxny;
                d(r,c)   = xy;
%                 if r==9 && c == 30 && i < 6 && j < 6
%                     figure,
%                     
%                     plot(sig1,'r');
%                     hold on
%                     plot(sig2,'b');
%                     title(['dot =',num2str(dotproduct1),' distance = ',num2str(d(r,c))])
%                 end
% %                 
                if ep~=0
                    d(r,c)=d(r,c)-ep;
                    d(r,c)=-d(r,c);
                    if d(r,c)>=0
                        d(r,c)=1;
                    else
                        d(r,c)=0;
                    end
                end
                if ep == 0
                    d(r,c)=abs(d(r,c));
                end
                if isnan(d(r,c)) == 0
                    p(r,c)=p(r,c)+d(r,c);
                end
                
            end
        end
    end
end

p=(p*2)/(Nsig*(Nsig-1));
p=p+1;
p=p.*double(mask);
% p=medfilt2(p);