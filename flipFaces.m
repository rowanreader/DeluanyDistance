function tri2 = flipFaces( tri, faces, vertices )
tri2 = tri;
% Change definition of two faces in the ConnectivityList
face1_verts = tri2.ConnectivityList(faces(1),:);
face2_verts = tri2.ConnectivityList(faces(2),:);
chgVerts = face1_verts(face1_verts~=vertices(2));
face1_verts(face1_verts==chgVerts(1)) = vertices(1);
face2_verts(face2_verts==chgVerts(2)) = vertices(2);
% Create updated triangulation object
ConnectivityList = tri2.ConnectivityList;
ConnectivityList(faces(1),:) = face1_verts;
ConnectivityList(faces(2),:) = face2_verts;
tri2 = triangulation(ConnectivityList,tri2.Points);
end