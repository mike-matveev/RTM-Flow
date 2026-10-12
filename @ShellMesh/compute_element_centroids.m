function obj = compute_element_centroids(obj)
%% compute_element_centroids  Mean 3D position of each element's nodes.
%
% Returns centroids as [num_elements x 3].  These 3D centroids are used
% for evaluating the permeability tensor K(centroid), exactly as the 2D
% code used 2D centroids — the user's K function simply accepts a 3-vector.

obj.centroids = zeros(obj.num_elements, 3);

for i = 1:obj.num_elements
    obj.centroids(i,:) = mean(obj.nodes(obj.elements(i,:), :), 1);
end

end
