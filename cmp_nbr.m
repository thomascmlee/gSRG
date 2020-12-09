function nbr_cells = cmp_nbr( seed_set, adj_mat, n, labeled_cells )
% Computes neighboring cells of a seed/growing region set.
%
% Args:
%   seed_set: Seed/growing region set consisting of a number of cells.
%   adj_mat: Adjacent matrix measuring distance between pairs of cells.
%   n: Total number of cells.
%   labeled_cells: Cells that have already been labeled.
%
% Returns:
%   nbr_cells: Neighboring cells of the seed/growing region set.

% Excludes seed_set itself and labeled_cells since they cannot be
% neighbors.
diff_set = setdiff(1:n, [seed_set labeled_cells]);

% Computes sum of each row.
nbr_cells = diff_set(sum(adj_mat(diff_set, seed_set), 2)>0);

end
