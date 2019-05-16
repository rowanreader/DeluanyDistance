% Insert Vertex into Triangulation
function [ tri2, is_new_face, is_new_vertex, new_vertex_ID ] = insert_vertex( tri, new_vertex, nearest_face, max_distance, is_new_face, is_new_vertex)
tri2 = tri;
% Adjacent triangles that have a dihedral angle that deviates from pi by an
% angle greater than filterangle are preserved (no Delaunay flipping is performed) 
% This conserves the previous triangulation shape. 
filterangle = pi/3;
% Do nothing if new_vertex is close to any existing vertex
nearest_face_def = tri2.ConnectivityList(nearest_face,:); % (1 x 3)   vertice IDs 
nearest_face_vertices = tri2.Points(nearest_face_def,:);  % (3 x xyz) vertice coordinates  
[distance,I] = min(sqrt(sum(bsxfun(@minus,new_vertex,nearest_face_vertices).^2,2)));
nearest_vertex_ID = nearest_face_def(I);
if distance<=max_distance
    new_vertex_ID = nearest_vertex_ID;
    return
end
% Three cases:
%   - Point is on a boundary edge
%   - Point is on an edge inside the mesh
%   - Point is inside a triangle (not on an edge)
% To create new faces, the old face is copied and vertice ID are replaced 
% to maintain the normal orientation of the face. 
% Change to matrix representation to be able to edit the triangulation
faces = tri2.ConnectivityList; % (#faces x 3)
vertices = tri2.Points;        % (#vertices x 3)
nVertices = size(tri2.Points,1);
% Distance to edge
[edgeDist,pt] = distance_to_edges(nearest_face_def,vertices,new_vertex,tri2.faceNormal);
if ~isnan(edgeDist) && abs(edgeDist) <= max_distance 
    % Point is on an edge
    bary = tri2.cartesianToBarycentric(nearest_face,pt); % (1 x 3)
    [~,I] = sort(bary); % ascending: first entry is the small one
    faceDef = faces(nearest_face,:); % (1 x 3)
    edge_verts = faceDef(I(2:3));    % (1 x 2) two vertices of edge
    edge_faces = tri2.edgeAttachments(edge_verts); % (1 x 1) cell array 
    edge_faces = edge_faces{1};      % (#edgeFaces x 1) faces at the edge
    
    FE = tri2.featureEdges(filterangle);
    if ismember(edge_verts,FE,'rows') || ismember(flip(edge_verts),FE,'rows')
        % Otherwise triangles with small angles will be produced
        new_vertex = pt;
    end
        
    if numel(edge_faces)==1 
        % Boundary edge
        
        % Add 1 new vertex
        vertices = [vertices;new_vertex];  % (#vertices x 3)
        is_new_vertex = [is_new_vertex;1]; % (#vertices x 1)
        new_vertex_ID = nVertices+1;       % scalar integer
        
        % Add 2 new faces
        new_faces = [faces(nearest_face,:); faces(nearest_face,:)]; % (2 x 3) copy old face 
        new_faces(1,I(2)) = new_vertex_ID; % (2 x 3) insert new vertex
        new_faces(2,I(3)) = new_vertex_ID; % (2 x 3) insert new vertex
        faces = [faces;new_faces];         % (#faces x 3)
        is_new_face = [is_new_face;1;1];   % (#faces x 1)
        
        % Remove 1 old face
        faces(nearest_face,:) = [];        % (#faces x 3)
        is_new_face(nearest_face) = [];    % (#faces x 1)
        
        % Indices for Delaunay check
        nFaces = size(faces,1);            % scalar integer
        delaunay_check = nFaces-1:nFaces;  % (2 x 1) vector
        
    elseif numel(edge_faces)==2 
    % Two attached faces on edge
        
        % Opposite vertices
        edge_faces_def = faces(edge_faces,:);      % (2 x 3)
        Lia = ismember(edge_faces_def,edge_verts); % logical (2 x 3), true where vertice IDs are on edge 
        replace_face1 = find(Lia(1,:));            % (1 x 2) position of edge's vertice IDs in edge face 1  
        replace_face2 = find(Lia(2,:));            % (1 x 2) position of edge's vertice IDs in edge face 2  
        
        % Add 1 new vertex
        vertices = [vertices;new_vertex];  % (#vertices x 3)
        new_vertex_ID = size(vertices,1);  % scalar integer
        is_new_vertex = [is_new_vertex;1]; % (#vertices x 1)
        
        % Add 4 new faces
        new_faces = [            % (4 x 3) copy old faces
            edge_faces_def(1,:)
            edge_faces_def(1,:)
            edge_faces_def(2,:)
            edge_faces_def(2,:)
            ];
        new_faces(1,replace_face1(1)) = new_vertex_ID; % (4 x 3) insert new vertex
        new_faces(2,replace_face1(2)) = new_vertex_ID; 
        new_faces(3,replace_face2(1)) = new_vertex_ID; 
        new_faces(4,replace_face2(2)) = new_vertex_ID; 
        
        faces = [faces;new_faces];           % (#faces x 3)
        is_new_face = [is_new_face;1;1;1;1]; % (#faces x 1)
        
        % Remove 2 old faces
        faces(edge_faces,:) = [];     % (#faces x 3)
        is_new_face(edge_faces) = []; % (#faces x 1)
        
        % Indices for Delaunay check
        nFaces = size(faces,1);           % scalar integer
        delaunay_check = nFaces-3:nFaces; % (1 x 4) vector
        
    else
        error('More than two faces (%s) on an edge. This case is not implemented yet.',mat2str(edge_faces))
    end    
else 
    % Point is not on an edge
    
    % Add 1 new vertex
    vertices = [vertices;new_vertex];  % (#vertices x 3)
    is_new_vertex = [is_new_vertex;1]; % (#vertices x 1)
    new_vertex_ID = nVertices+1;       % scalar integer
    
    % Add 3 new faces
    nearest_face_def = faces(nearest_face,:); % (1 x 3)
    new_faces = repmat(nearest_face_def,3,1); % (3 x 3) copy old face
    new_faces(1,1) = new_vertex_ID;           % insert new vertex ID
    new_faces(2,2) = new_vertex_ID;
    new_faces(3,3) = new_vertex_ID;
    faces = [faces;new_faces];         % (#faces x 3)
    is_new_face = [is_new_face;1;1;1]; % (#faces x 1)
    
    % Remove 1 old face
    faces(nearest_face,:) = [];     % (#faces x 1)
    is_new_face(nearest_face) = []; % (#faces x 1)
    
    % Indices for Delaunay check
    nFaces = size(faces,1);           % scalar
    delaunay_check = nFaces-2:nFaces; % (1 x 3) vector
end
% Change to triangulation class representation again
tri2 = triangulation(faces,vertices);
% Restore Delaunay conditions
for k = delaunay_check
    [tri2, changed_faces] = makeDelaunay(tri2, k, filterangle);  
    is_new_face = is_new_face | changed_faces; % (#faces x 1)
end
end