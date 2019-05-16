function [D,P,F,V] = distance_to_vertices(faces,vertices,qPoint,normals)
% find nearest vertex
[D,nearestVertexID] = min(sum(bsxfun(@minus,vertices,qPoint).^2,2),[],1);
D = sqrt(D);
P = vertices(nearestVertexID,:); % (1 x 3)
V = nearestVertexID;
% find faces that belong to the vertex
connectedFaces = find(any(faces==nearestVertexID,2)); % (#connectedFaces x 1) face indices
assert(length(connectedFaces)>=1,'Vertex %u is not connected to any face.',nearestVertexID)
F = connectedFaces(1);
n = normals(connectedFaces,:); % (#connectedFaces x 3) normal vectors
% scalar product between distance vector and normal vectors
coefficients = dot2(n,qPoint-P);
sgn = signOfLargest(coefficients);
D = D*sgn;
end