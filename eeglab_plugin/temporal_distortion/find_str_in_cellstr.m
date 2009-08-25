function [idx] = find_str_in_cellstr( cellstr, mystr );
  idx=0;
  for i=1:size(cellstr,2)
    if strcmp(mystr, char(cellstr{i}))
      idx=i;
      break;
    end;
  end;
  
  return;