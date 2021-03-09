function mat_out = bush_recurse(mat_in, width_in)
  mat_out = mat_in;
  row_last = mat_in(size(mat_in, 1) , :);
  i = 2;
  while i <= width_in
    if (row_last(i-1) > (row_last(i) + 1)  )
      row_new = row_last;
      row_new(i) = row_new(i) + 1;
      row_new(i-1) = row_new(i-1) - 1;
      mat_out = [mat_out; row_new];
      row_last = row_new;
      i = 2;
    elseif (i + 1 <= width_in)
        if ((row_last(i-1) > row_last(i) ) && (row_last(i) > row_last(i+1) ) ) 
            row_new = row_last;
            row_new(i+1) = row_new(i+1) + 1;
            row_new(i-1) = row_new(i-1) - 1;
            mat_out = [mat_out; row_new];
            row_last = row_new;
            i = 2;
        else
            i = i + 1;
        end
    else
        i = i + 1;
    end
        
  end
  

  
end
