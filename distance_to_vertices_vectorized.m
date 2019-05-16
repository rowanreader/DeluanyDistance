% Vectorized Distance Functions 
%  (can process more than one point on the same mesh)
function [D,P,F,V] = distance_to_vertices_vectorized(faces,vertices,points,normals)
% Requires Statistics Toolbox 
% [D,I] = pdist2(vertices,points, 'euclidean', 'Smallest',1); % (1 x #points)
% D = D'; % (#points x 1)
% I = I'; % (#points x 1)
D = sum(bsxfun(@minus,permute(vertices,[3,2,1]),points).^2,2); % (#points x 1 x #vertices)
[D,I] = min(D,[],3); % (#points x 1)
D = sqrt(D);         % (#points x 1)
P = vertices(I,:);   % (#points x 3)
V = I;
% find faces that belong to the vertex
inds = I; % (#points x 1)
inds = permute(inds,[3,2,1]); % (1 x 1 x #points)
inds = bsxfun(@eq,faces,inds); % (#faces x 3 x #points)
inds = any(inds,2); % (#faces x 1 x #points) logical indices which faces belong to the nearest edge of a query point
inds = permute(inds,[1,3,2]); % (#faces x #points)
inds = num2cell(inds,1); % (1 x #points) cell array with (#faces x 1) logical indices
n = cellfun(@(x) normals(x,:), inds, 'UniformOutput',false)'; % (#points x 1) cell array with (#connectedFaces x 3) normal vectors
F = cellfun(@(x) find(x,1), inds)'; % (#points x 1)
% scalar product between distance vector and normal vectors
coefficients = cellfun(@dot2, n, num2cell(points-P,2), 'UniformOutput',false);
sgn = cellfun(@signOfLargest,coefficients);
D = D.*sgn;
end