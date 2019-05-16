function [tri2,highlight_faces2,faces_2To1,vertices_2To1] = subSurface( tri, vertex, highlight_faces, iterations )
% enlarge surface by vertex attachments
nFaces = size(tri.ConnectivityList,1);
nVertices = size(tri.Points,1);
faceL    = ind2logical([],nFaces);
verticeL = ind2logical([],nVertices);
verticeL(vertex) = true;
for x = 1:iterations
    faceIDs = tri.vertexAttachments(find(verticeL))'; %#ok<FNDSB> no array indexing 
    %faceIDs = cell2mat(faceIDs);
    faceIDs = [faceIDs{:}];
    faceL(faceIDs) = true;
    vertexIDs = tri.ConnectivityList(faceL,:);
    verticeL(vertexIDs) = true;
end
% faceL = any(tri.ConnectivityList==vertex,2);
% for x = 1:iterations-1
%     verts1 = tri.ConnectivityList(faceL,:);
%     faceL = any(ismember(tri.ConnectivityList,verts1),2);
% end
[connectedVertices,~,newIDs] = unique(tri.ConnectivityList(faceL,:));
faces2 = reshape(newIDs,[],3);               % (#connectedFaces    x 3) new face connectivity list
vertices2 = tri.Points(connectedVertices,:); % (#connectedVertices x 3) new vertice list
tri2 = triangulation(faces2,vertices2);
highlight_faces2 = ind2logical(highlight_faces,nFaces);
highlight_faces2 = highlight_faces2(faceL);
highlight_faces2 = find(highlight_faces2);
faces_2To1 = find(faceL);
vertices_2To1 = connectedVertices;
end
