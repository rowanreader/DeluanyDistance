function [D,P,F] = distance_to_surfaces_vectorized(faces,vertices,points,normals)
r1 = vertices(faces(:,1),:);   % (#faces x 3) % 1st vertex of every face 
r2 = vertices(faces(:,2),:);   % (#faces x 3) % 2nd vertex of every face 
r3 = vertices(faces(:,3),:);   % (#faces x 3) % 3rd vertex of every face 
qp = permute(points,[3,2,1]);  % (1 x 3 x #points) query points
vq = bsxfun(@minus,qp,r1);     % (#faces x 3 x #points) 
D = dot2(vq,normals);          % (#faces x 1 x #points) distance to surface
rD = bsxfun(@times,normals,D); % (#faces x 3 x #points) vector from surface to query point 
P = bsxfun(@minus,qp,rD);      % (#faces x 3 x #points) nearest point on surface; can be outside triangle 
% find barycentric coordinates (query point as linear combination of two edges) 
r31r31 = sum((r3-r1).^2,2);    % (#faces x 1)
r21r21 = sum((r2-r1).^2,2);    % (#faces x 1)
r21r31 = dot(r2-r1,r3-r1,2);   % (#faces x 1)
r31vq = dot2(r3-r1,vq);        % (#faces x 1 x #points)
r21vq = dot2(r2-r1,vq);        % (#faces x 1 x #points)
d = r31r31.*r21r21 - r21r31.^2; % (#faces x 1)
bary = NaN(size(faces,1), 2, size(points,1)); % (#faces x 3 x #points) 
bary(:,1,:) = bsxfun(@rdivide, bsxfun(@times,r21r21,r31vq) - bsxfun(@times,r21r31,r21vq), d); 
bary(:,2,:) = bsxfun(@rdivide, bsxfun(@times,r31r31,r21vq) - bsxfun(@times,r21r31,r31vq), d); 
bary(:,3,:) = 1 - bary(:,1,:) - bary(:,2,:);  % (#faces x 3 x #points) 
% exclude intersections that are outside the triangle
D( any(bary<=0,2) | any(bary>=1,2) ) = NaN;  % (#faces x 1 x #points)
D( abs(d)<=eps, :, : ) = NaN;
% find nearest face for every query point
[~,I] = min(abs(D),[],1); % (1 x 1 x #points)
I = squeeze(I); % (#points x 1)
D = D(sub2ind(size(D),I,ones(length(I),1),(1:length(I))'));
D = squeeze(D); % (#points x 1)
P = permute(P,[2,1,3]); % (3 x #faces x #points)
sz = [size(P) 1 1];
P = P(:,sub2ind(sz(2:3),I,(1:length(I))')); % (3 x #points)
P = P'; % (#points x 3)
F = I;  % (#points x 1)
end