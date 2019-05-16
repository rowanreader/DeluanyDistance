function [D,P,F] = distance_to_edges(faces,vertices,qPoint,normals)
% Point-point representation of all edges
edges = [faces(:,[1,2]); faces(:,[1,3]); faces(:,[2,3])]; % (#edges x 2) vertice IDs
% Intersection between tangent of edge lines and query point
r1 = vertices(edges(:,1),:);   % (#edges x 3) first point of every edge
r2 = vertices(edges(:,2),:);   % (#edges x 3) second point of every edge
t = dot( bsxfun(@minus,qPoint,r1), r2-r1, 2) ./ sum((r2-r1).^2,2); % (#edges x 1) location of intersection relative to r1 and r2
t(t<=0) = NaN; % exclude intersections not between the two vertices r1 and r2
t(t>=1) = NaN;
% Distance between intersection and query point
P = r1 + bsxfun(@times,(r2-r1),t); % (#edges x 3) intersection points
D = bsxfun(@minus,qPoint,P); % (#edges x 3)
D = sqrt(sum(D.^2,2));       % (#edges x 1)
[D,I] = min(D,[],1);         % (1 x 1)
P = P(I,:);
% find faces that belong to the edge
inds = edges(I,:);  % (1 x 2)
inds = permute(inds,[3,1,2]);  % (1 x 1 x 2)
inds = bsxfun(@eq,faces,inds); % (#faces x 3 x 2)
inds = any(inds,3);    % (#faces x 3)
inds = sum(inds,2)==2; % (#faces x 1) logical indices which faces belong to the nearest edge of the query point
F = find(inds,1);
n = normals(inds,:); % (#connectedFaces x 3) normal vectors
% scalar product between distance vector and normal vectors
coefficients = dot2(n,qPoint-P); % (#connectedFaces x 1)
sgn = signOfLargest(coefficients);
D = D*sgn;
end
