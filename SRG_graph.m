function [init_sets, labeled_cells] = SRG_graph( init_sets, cell_log_intensity, cell_area, n, adj_mat, invalid, verbose, print_n )
% Graph-based seeded region growing.
%
% Args:
%   init_sets: (cell) Initial seed/growing region sets represented by indices
%     of Voronoi cells.
%   cell_log_intensity: Log intensity of Voronoi cells.
%   cell_area: Area of Voronoi cells.
%   n: Total number of Voronoi cells.
%   adj_mat: Adjacent matrix that gives the connectivity of Voronoi cells.
%   invalid: Indices of invalid Voronoi cells; invalid is due to that the
%     area is zero or the Voronoi cell is on the boundary of the domain.
%   verbose: Whether print out the status.
%   print_n: Print every n iterations.
%
% Returns:
%   init_sets: (cell) The region sets after finish growing.
%   labeled_cells: The Voronoi cells that have been assigned to one of the
%     region sets.

if nargin==6
    verbose = false;
end

% Initialize seed/growing region sets.
m = length(init_sets);
region_log_intensity = zeros(1, m);
region_area = zeros(1, m);
labeled_cells = [];

% Neighboring Voronoi cells of each growing region set.
nbr_cells = cell(m, 1);
for j = 1:m
    % Note that init_sets{j} is a row vector.
    % Update labeled_cells.
    labeled_cells = [labeled_cells init_sets{j}];
    region_area(j) = sum(cell_area(init_sets{j}));
    region_log_intensity(j) = log(sum(cell_area(init_sets{j}).*...
        exp(cell_log_intensity(init_sets{j})))/region_area(j));
end

% Compute sets of neighbors for each growing region set.
for j = 1:m
    nbr_cells{j} = cmp_nbr( init_sets{j}, adj_mat, n, [labeled_cells invalid] );
end

% Each row of min_c_pairs contains two elements: the min diff between the 
% growing region and its neighboring cells, and the index of the cell with
% min diff.
% Keeping a record of them significantly speeds up the computation.
min_c_pairs = zeros(m, 2);
for j = 1:m
    % If region j has neighbors.
    if ~isempty(nbr_cells{j})
        [min_c, min_index] = min(abs(cell_log_intensity(nbr_cells{j})-region_log_intensity(j)));
        min_c_pairs(j, 1) = min_c;
        min_c_pairs(j, 2) = nbr_cells{j}(min_index);
    else
        min_c_pairs(j, 1) = realmax;
        min_c_pairs(j, 2) = 0;
    end
end

% Get the number of Voronoi cells that have been assigned to the initial 
% seed/growing region sets.
n_init = length(labeled_cells);
% Get the number of invalid Voronoi cells.
n_invalid = length(invalid);
% Get the number of remaining Voronoi cells.
n_remain = n-n_init-n_invalid;

% Start the iteration to assign the remaining Voronoi cells. 
for i = 1:n_remain
    
    if verbose
        if mod(i, print_n)==0
            disp(['Iteration ', num2str(i), '...'])
        end
    end
    
    % Select the pair of a growing region set and a neighboring Voronoi 
    % cell with the smallest value of the criterion, where j refers to 
    % regions, and k refers to cells.
    [~, min_j] = min(min_c_pairs(:, 1));
    min_k = min_c_pairs(min_j, 2);
    if min_k == 0
        disp('Warning: There are isolated cells so the algorithm terminates prematurely.')
        disp(['There are in total ', num2str(n_remain - i + 1), ' isolated cells.'])
        break
    end
    
    % Update the selected growing region set.
    labeled_cells = [labeled_cells min_k];
    init_sets{min_j} = [init_sets{min_j} min_k];
    region_log_intensity(min_j) = log((exp(region_log_intensity(min_j))*...
        region_area(min_j)+exp(cell_log_intensity(min_k))*cell_area(min_k))/...
        (region_area(min_j)+cell_area(min_k)));
    region_area(min_j) = region_area(min_j)+cell_area(min_k);
    % Add the set of neighbors for cell min_k to that of region min_j, and
    % make sure the uniqueness of the elements.
    nbr_cells{min_j} = unique([nbr_cells{min_j} cmp_nbr( min_k, adj_mat, n, [labeled_cells invalid] )]);
   
    for j = 1:m
        if ismember(min_k, nbr_cells{j})
            % Delete cell min_k from the set of neighbors for each growing 
            % region set.
            nbr_cells{j} = setdiff(nbr_cells{j}, min_k);
            if min_c_pairs(j, 2)==min_k
                % Update min_c_pairs for region j.
                if ~isempty(nbr_cells{j})
                    [min_c, min_index] = min(abs(cell_log_intensity(nbr_cells{j})-region_log_intensity(j)));
                    min_c_pairs(j, 1) = min_c;
                    min_c_pairs(j, 2) = nbr_cells{j}(min_index);
                else
                    min_c_pairs(j, 1) = realmax;
                    min_c_pairs(j, 2) = 0;
                end
            end
        end
    end
    
end
        
end
