function [D,P,F] = distance_to_surfaces(faces,vertices,point,normals)
r1 = vertices(faces(:,1),:);   % (#faces x 3) % 1st vertex of every face
r2 = vertices(faces(:,2),:);   % (#faces x 3) % 2nd vertex of every face
r3 = vertices(faces(:,3),:);   % (#faces x 3) % 3rd vertex of every face
vq = bsxfun(@minus,point,r1);  % (#faces x 3)
D = dot(vq,normals,2);         % (#faces x 1) distance to surface
rD = bsxfun(@times,normals,D); % (#faces x 3) vector from surface to query point
P = bsxfun(@minus,point,rD);   % (#faces x 3) nearest point on surface; can be outside triangle
% find barycentric coordinates (query point as linear combination of two edges)
r31r31 = sum((r3-r1).^2,2);    % (#faces x 1)
r21r21 = sum((r2-r1).^2,2);    % (#faces x 1)
r21r31 = dot(r2-r1,r3-r1,2);   % (#faces x 1)
r31vq = dot(r3-r1,vq,2);       % (#faces x 1)
r21vq = dot(r2-r1,vq,2);       % (#faces x 1)
d = r31r31.*r21r21 - r21r31.^2;               % (#faces x 1)
bary = NaN(size(faces,1),3);                  % (#faces x 3)
bary(:,1) = (r21r21.*r31vq-r21r31.*r21vq)./d; % (#faces x 3)
bary(:,2) = (r31r31.*r21vq-r21r31.*r31vq)./d; % (#faces x 3)
bary(:,3) = 1 - bary(:,1) - bary(:,2);        % (#faces x 3)
% tri = triangulation(faces,vertices);
% bary = tri.cartesianToBarycentric((1:size(faces,1))',P); % (#faces x 3)
% exclude intersections that are outside the triangle
D( abs(d)<=eps | any(bary<=0,2) | any(bary>=1,2) ) = NaN;  % (#faces x 1)
% find nearest face for query point
[~,I] = min(abs(D),[],1); % (1 x 1)
D = D(I);       % (1 x 1)
P = P(I,:);     % (1 x 3)
F = I;          % (1 x 1)
end
