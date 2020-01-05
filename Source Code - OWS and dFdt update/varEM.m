
function [x]=varEM(fstr,varname,numswitch)
% function for finding specified variable from EM config file (produced automatically using MATLAB)

%% Find Variable
ind=strfind(fstr,[varname,'=']);

%% Loop incase of more than 1 start found (i.e. t=) 
if numel(ind) > 1
   for j=1:numel(ind)
        k=strcmp(fstr(ind(j)-3),'!');
        if k == 1
            indi=ind(j);
        end
   end
ind=[];
ind=indi;
end

%% Find numerical value (as string)
first_ind=ind+numel(varname)+1;
for i =1:[numel(fstr)-first_ind]
    k=strcmp(fstr(first_ind+i),'!');
    if k==1
    last_ind=first_ind+i-1;
    break
    end
end
x = fstr(first_ind:last_ind);

%% Switch if variable is value of popupbox/checkbox (i.e. not a string)
if numswitch == 1
    x=str2num(x);
end