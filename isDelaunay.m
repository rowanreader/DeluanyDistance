function isDelaunay = isDelaunay( tri, faces, vertices )
% Checks if two faces (that share an edge) meet the Delaunay condition
%    - tri:      triangulation object  
%    - faces:    (2 x 1) face IDs
%    - vertices: (2 x 1) vertice IDs of opposite vertices
[CC,r] = tri.circumcenter(faces); % CC: (2 x xyz) centers;  r: (2 x 1) radius
distance = sqrt(sum((tri.Points(vertices,:)-CC).^2,2)); % (2 x 1) distances
isDelaunay = all( distance > r+eps );
end
