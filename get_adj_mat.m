function adj_mat = get_adj_mat( E, n )
% Gets adjacent matrix.
% 
% Args:
%   E: A matrix each row of which represents the two vertices of an edge.
%   n: Total number of cells.
%
% Returns:
%   adj_mat: The adjacent matrix, where 1 represents that there is an edge,
%   and 0 represents that there is no edge.

adj_mat = false(n);
for i = 1:size(E, 1)
    adj_mat(E(i, 1), E(i, 2)) = true;
    adj_mat(E(i, 2), E(i, 1)) = true;
end

end
