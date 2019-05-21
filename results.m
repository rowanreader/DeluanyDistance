% calculates surface points, mainly

function [dist, spts, f2,v2,vID, fID] = results(file)

addpath("C:\Users\Jacqueline\Documents\MATLAB\Jacqueline\Breast Demarkation");
[v,f,n,~] = stlRead('S1.stl');
pts = inputPts(file); % output or output2
[m,~] = size(pts);
[dist, spts, f2, v2, vID, fID] = point2trimesh('Faces', f, 'Vertices', v, 'QueryPoints',pts, 'MaxDistance', 1);

faceNorms = zeros(m,3);
id = zeros(m,1);
[o,~] = size(n);

% find which face each surface point belongs to
% use dot product of (spts - v, n)

for i = 1:m
    minN = 1;
    for j = 1:o
        temp = dot(spts(i,:)-v(f(j,1),:),n(j,:)); % check if on plane
        if abs(temp) < minN && inTriangle(spts(i,:), v(f(j,1),:), v(f(j,2),:), v(f(j,3),:)) % check if better and in face
            %tt = inTriangle(spts(i,:), v(f(j,1),:), v(f(j,2),:), v(f(j,3),:));
            id(i) = j;
            faceNorms(i,:) = n(j,:);
        end
    end
end
dOut = fopen('distances.txt','w');
for i=1:m
    fprintf(dOut, '%f\n', dist(i));
end
fclose(dOut);

surfPts = fopen('surface points.txt','w');
for i=1:m
    fprintf(surfPts, '%f %f %f\n', spts(i,:));
end
fclose(surfPts);

nOut = fopen('normals.txt','w');
for i=1:m
    fprintf(nOut, '%d %f %f %f\n', id(i), faceNorms(i,1), faceNorms(i,2), faceNorms(i,3));
end
fclose(nOut);
