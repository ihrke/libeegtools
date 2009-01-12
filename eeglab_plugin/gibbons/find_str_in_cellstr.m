function [idx] = find_str_in_cellstr( cellstr, mystr );
  idx=0;
  for i=1:size(cellstr,1)
    if strcmp(mystr, cellstr{i})
      idx=i;
      break;
    end;
  end;
  
  return;