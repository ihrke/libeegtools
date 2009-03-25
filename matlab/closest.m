function [ index ] = closest( v, v2 )
%closest return index of closest entry in vector v to scalar x (abs)
index=[];

if size(v2,1)>size(v2,2)
    v2 = v2';
end;
for x=v2
    [tmp idx] = min(abs(v-x));
    index = [index idx];
end;
