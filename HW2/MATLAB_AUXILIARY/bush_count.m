function count = bush_count(tree_count, r, bush)

  part_mat = partition_matrix(r, bush); %% should yield all possible sub-tree configurations
  count = 0;

  for i = 1:size(part_mat, 1)
    local_row = part_mat(i, :);
    [tree_occur, tree_id] = groupcounts(local_row');
    tree_id = tree_id + 1;
    acc = 1;
    for j = 1:length(tree_id)
        n = tree_count(tree_id(j)) + tree_occur(j) - 1;
        acc =acc*( nchoosek( n, tree_occur(j) ) );
    end
    count = count + acc;
  end

end
