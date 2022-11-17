function o=orthog(x,y)

% RH 1999

% orthog x wrt y

if(size(x,1)==1)
x=x';
end
if(size(y,1)==1)
y=y';
end

o=x-y*pinv(y)*x;