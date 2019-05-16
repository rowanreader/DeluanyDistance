% Restore Delaunay Conditions
function [tri, changed_faces] = makeDelaunay( tri, face, filterangle )
% Ensures that the Delaunay condition around a given face is met.
% If an edge has to be flipped, the resulting changed faces are checked, too.
% TODO Currently, more Delaunay checks than necessary are performed
changed_faces   = false(size(tri.ConnectivityList,1),1); % (#faces x 1) logical vector which faces that have been flipped
unchecked_faces = false(size(tri.ConnectivityList,1),1); % (#faces x 1) logical vector which faces still have to be checked
unchecked_faces(face) = true;
while any(unchecked_faces)
    current_face = find(unchecked_faces,1);
    unchecked_faces(current_face) = false;
    [tri, changed_faces2] = makeDelaunay_neighbors( tri, current_face, filterangle );
    unchecked_faces = unchecked_faces | changed_faces2;
    changed_faces   = changed_faces   | changed_faces2;
end
end
