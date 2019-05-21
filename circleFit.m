
% pulls out points that are closer to being a sphere instead of chest wall
% fits circle to it
% gets points along circle surface
function unrotated = circleFit
global dirNum;
t = 10;
plot = 1;
spts = inputPts('surface points.txt');
if plot == 1
    figure;
    plot3(spts(:,1), spts(:,2), spts(:,3),'b.');
    axis equal;
end
% only get the start of each line?
[n,~] = size(spts);
dirCount = 1;
index = 1;
thresh = 6; % first thresh points THIS IS WHAT YOU CHANGE
circPts = zeros(36*thresh,3); % num of lines * num points
for i=1:n
    if dirCount < thresh
        circPts(index,:) = spts(i,:);
        index = index + 1;
    end
    if dirCount == dirNum
        dirCount = 0;
    end
    dirCount = dirCount + 1;
end
if plot == 1
    figure;
    plot3(circPts(:,1), circPts(:,2), circPts(:,3),'r.');
    %axis equal;
    hold on;
end



% need to fit a sphere to circPts
[center, rad] = sphereFit(circPts);
rad = rad - 25;
[x,y,z] = sphere();
sphX = x*rad + center(1);
sphY = y*rad + center(2);
sphZ = z*rad + center(3);
if plot == 1
    surf(sphX, sphY, sphZ);
    axis equal;
    alpha 0.3;
    hold on;
end



pts = inputPts('sphere_points_leftBreast_S1.txt');
[~,circRadius, circCenter, normal] = vectors(pts); % literally just to get the radius and center of the inner circle
normal = normal/(norm(normal));
% rotate about x, then y => also make 'unrotate' matrices
a = normal(1);
b = normal(2);
c = normal(3);
dir = sqrt(normal(2).^2 + normal(3).^2);

v = sqrt(b^2 + c^2);
cos_t = c/v;
sin_t = b/v;

rotate_xi = [1,0,0;...
    0,cos_t, sin_t;...
    0, -sin_t, cos_t];

rotate_yi = [dir, 0, a;...
    0,1,0;...
    -a, 0, dir];


% get points on circle to start with (10 degrees apart) - literally just
% copy from vectors.m
num = 360/t;
% we have center, normal, radius
%dirNum = 30; % how many dots/line
space = 5; % separation of dots
temp = zeros(dirNum,3); % holds points for a given line
lines = zeros(num*dirNum,3); % holds points for all lines
for i = 0:num
    degree = i*t; % every t degrees
    % find point on circle at degree
    x = circRadius*sind(degree);
    y = circRadius*cosd(degree);    
    % find perpindicular to tangent (literally just treat as direction
    % vector)
    vector = [x,y];
    vector = vector/norm(vector); % normalize
    % find dirNum points along line, curving around circle
    % every 10 degrees out
    % find point on line
    % find vector from point to center
    % find point sphere's radius away along vector
    for j = 1:dirNum
        vecDist = space * j + circRadius; % lenght of vector to point
        vecPoint = vector*vecDist; % actual point;
        tempDirVec = [vecPoint,0] - center; % direction vector from point to center, z is 0, will rotate later
        tempDirVec = tempDirVec/norm(tempDirVec);
        tempPoint = center + (tempDirVec * rad); % actual point
        % since at origin and in line with z axis, all z coordinates are 0
        temp(j,:) = tempPoint;
    end
    lines(i*dirNum+1:(i+1)*dirNum,:) = temp; 
end

unrotated =  (circCenter' + (-1*rotate_yi*(rotate_xi*lines'))); 
% -1 because you screwed something up with the flips
%unrotated =  circCenter' + ((rotate_xi*lines'));

unrotated = unrotated';
% unrotated(1) = unrotated(1) + 20;
% unrotated(2) = unrotated(2) + 10;
% unrotated(3) = unrotated(3) + 10;
if plot == 1    
    plot3(unrotated(:,1), unrotated(:,2), unrotated(:,3),'g.');
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
end


% write unrotated to file
file = fopen('output2.txt','w');

for i =1:(num*dirNum)
    fprintf(file,'%f %f %f\n', unrotated(i,1),unrotated(i,2),unrotated(i,3));
end
fclose(file);




