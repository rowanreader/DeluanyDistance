function [D,P,F] = distance_to_edges_vectorized(faces,vertices,points,normals)
euclid = @(A,dim) sqrt(sum(A.^2,dim));
% Point-point representation of all edges
edges = [faces(:,[1,2]); faces(:,[1,3]); faces(:,[2,3])]; % (#edges x 2) vertice IDs 
% Intersection between tangent of edge lines and query points
r1 = vertices(edges(:,1),:);   % (#edges x 3) first point of every edge 
r2 = vertices(edges(:,2),:);   % (#edges x 3) second point of every edge
qp = permute(points,[3,2,1]);  % (1 x 3 x #points) query points
t = bsxfun(@rdivide,...
        dot2(bsxfun(@minus,qp,r1), r2-r1),...
        sum((r2-r1).^2,2)); % (#edges x 1 x #points) location of intersection relative to r1 and r2 
t(t<=0) = NaN; % exclude intersections not between the two vertices r1 and r2  
t(t>=1) = NaN;
% Distance between intersection and query points
P = bsxfun(@plus,...
        r1,...
        bsxfun(@times,...
            (r2-r1),...
            t)); % (#edges x 3 x #points) intersection points
D = bsxfun(@minus,qp,P); % (#edges x 3 x #points) 
D = euclid(D,2);         % (#edges x 1 x #points) 
[D,I] = min(D,[],1);     % (1 x 1 x #points) 
D = squeeze(D);          % (#points x 1)
I = squeeze(I);          % (#points x 1)
P = permute(P,[2,1,3]);  % (3 x #edges x #points)
sz = [size(P) 1 1];
P = P(:,sub2ind(sz(2:3),I,(1:length(I))')); % (3 x #points)
P = P';                  % (#points x 3)
% find faces that belong to the edge
inds = edges(I,:);  % (#points x 2)
inds = permute(inds,[4,3,1,2]); % (1 x 1 x #points x 2)
inds = bsxfun(@eq,faces,inds);  % (#faces x 3 x #points x 2)
inds = any(inds,4);    % (#faces x 3 x #points)
inds = sum(inds,2)==2; % (#faces x 1 x #points) logical indices which faces belong to the nearest edge of a query point
inds = permute(inds,[1,3,2]); % (#faces x #points)
inds = num2cell(inds,1);      % (1 x #points) cell array with (#faces x 1) logical indices
n = cellfun(@(x) normals(x,:), inds, 'UniformOutput',false)'; % (#points x 1) cell array with (#connectedFaces x 3) normal vectors
F = cellfun(@(x) find(x,1), inds)';  % (#points x 1)
% scalar product between distance vector and normal vectors
coefficients = cellfun(@dot2, n, num2cell(points-P,2), 'UniformOutput',false);
sgn = cellfun(@signOfLargest,coefficients);
D = D.*sgn;
end