function [ r, ind ] = strisin( A, s )
%STRISIN Check if string s is in A (cell array), and return corr. bool.
%  * if ind is a cell array, return a bool for each element in s along with
%    the indices, where they have been found in A
%  * 
r = 0;
if(iscell(s))
    for el = 1:length(s)
        ind = strfind(A, s{el});
        for i = 1:length(ind)
            if(~isempty(ind{i})) 
                r = 1;
                break;
            end;
        end;
        
    end;
else % it is a string
    r=0;
    ind = strfind(A, s);
    i = 1;
    for i = 1:length(ind)
        if(~isempty(ind{i})) 
            r = 1;
            return;
        end;
    end;

end;

return;