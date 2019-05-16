function [tri, changed_faces] = makeDelaunay_neighbors( tri, face, filterangle )
% Ensures that the Delaunay condition between the given face and its neighbors is met
% After one flipping is performed, the corresponding faces are marked 
% unchecked and the function returns. 
% For every neighbor face, check Delaunay condition and flip edge if necessary
changed_faces = false(size(tri.ConnectivityList,1),1); % (#faces x 1) logical vector which faces that have been changed 
% Find neighbor faces and their vertices
attached_faces = tri.neighbors(face)';      % (3 x 1) neighbor face IDs
attached_faces(isnan(attached_faces)) = []; % (#attached_faces x 1) face IDs
attached_vertice_IDs = tri.ConnectivityList(attached_faces,:); % (#attached_faces x 3) vertice IDs
own_vertice_IDs = tri.ConnectivityList(face,:);       % (1 x 3) vertice IDs
lia = ismember(attached_vertice_IDs,own_vertice_IDs); % (#attached_faces x 3) logical index
attached_vertice_IDs_T = attached_vertice_IDs';       % (3 x #attached_faces) vertice IDs
opposite_vertice_IDs = attached_vertice_IDs_T(~lia'); % (#attached_faces x 1) vertice IDs
for neighborID = 1:size(attached_faces,1)
    face1 = face;
    face2 = attached_faces(neighborID);
    faces = [face1;face2]; % (2 x 1) face IDs
    vertex1 = opposite_vertice_IDs(neighborID);
    vertex2 = own_vertice_IDs(~ismember(own_vertice_IDs,attached_vertice_IDs(neighborID,:)));
    vertices = [vertex1;vertex2]; % (2 x 1) vertice IDs
    assert(length(vertices)==2)
    if ~isDelaunay( tri, faces, vertices )
        % Where sharp/acute/small angles occur, there is no Delaunay possible
        % TODO if points to insert are outside the surface, the feature
        % edges should be detected before any change of the mesh is performed
        % (This is not the case here, since we use only surface projection points) 
        FE = tri.featureEdges(filterangle);
        edge = intersect(tri.ConnectivityList(faces(1),:),tri.ConnectivityList(faces(2),:));
        sharpEdge = ismember(edge,FE,'rows') || ismember(flip(edge),FE,'rows');
        % don't flip if other neighbors have more than one common vertex to flip partner 
        neighbors = tri.neighbors(faces(1));
        neighbors(isnan(neighbors)) = [];
        neighbors(neighbors==faces(2)) = []; % neighbors of faces(1), without flip partner faces(2) 
        nNeighbors = arrayfun(@(neighbor) length(intersect(tri.ConnectivityList(neighbor,:),tri.ConnectivityList(faces(2),:))), neighbors);
        assert(~any(nNeighbors==0))
        commonVertices = any(nNeighbors>1);
        if ~sharpEdge && ~commonVertices
            tri2 = flipFaces( tri, faces, vertices );
            % Check if minimum angle increased by the flip
            if minimumAngle(tri2,faces)>minimumAngle(tri,faces) 
                tri = tri2;
                changed_faces(faces) = true;
                return
            end
        end
    end
end
end
