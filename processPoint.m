% Non-vectorized Distance Functions
%  (can process only one point)
function [D,P,F] = processPoint(faces,vertices,point,normals, distance_to_vertices,distance_to_edges,distance_to_surfaces, useSubSurface)
d = NaN(3,1); % (distanceTypes x 1)
p = NaN(3,3); % (distanceTypes x xyz)
f = NaN(3,1); % (distanceTypes x 1)
% find nearest vertice
[d(1),p(1,:),f(1),v] = distance_to_vertices(faces,vertices,point,normals);
% d:  (1 x 1) signed distance to surface
% p:  (1 x 3) corresponding point on surface
% v:  (1 x 1) nearest vertex
% connectedFaces: (#connectedFaces x 1) face indices
if useSubSurface
    [tri2,~,faces_2To1] = subSurface( triangulation(faces,vertices), v, [], 2 );
    faces2 = tri2.ConnectivityList;
    vertices2 = tri2.Points;
    normals = normals(faces_2To1,:);
else
    faces2   = faces;
    vertices2 = vertices;
end
% find nearest point on all edges
[d(2),p(2,:),f(2)] = distance_to_edges(faces2,vertices2,point,normals);
% find nearest point on all surfaces
[d(3),p(3,:),f(3)] = distance_to_surfaces(faces2,vertices2,point,normals);
if useSubSurface
    % translate back f(2) and f(3)
    f(2:3) = faces_2To1(f(2:3));
end
% find minimum distance type
[~,I] = min(abs(d),[],1);
D = d(I);
P = p(I,:);
F = f(I);
end
