function ang = angles
global dirNum;
file = fopen('normals.txt','r');

format long;
data = fscanf(file, '%f');
[n,~] = size(data);
m = n/4;
normal = zeros(m,4);
count = 1;
% should be in order, starting with 1st line at angle 0, all dirNum points 
% going out,  going clockwise 
for i = 1:m
    normal(i,1) = data(count);
    count = count + 1;
    normal(i,2) = data(count);
    count = count + 1;
    normal(i,3) = data(count);
    count = count + 1;
    normal(i,4) = data(count);
    count = count + 1;
end
fclose(file);
% for each set of dirNum points (1 line), compare angle between current and
% next
% will only go up to dirNum -1
ang = zeros(m,1);
count = 1;
angleOut = fopen('angles.txt','w');
% first value in normal is index of face - disregard
for j = 1:m
    if count == dirNum
        count = 1;
        continue;
    else
        a = normal(j,2:4);
        b = normal(j+1,2:4);
        if norm(a) ~= 0 && norm(b) ~= 0            
            ang(j) = acosd(dot(a,b)/(norm(a)*norm(b)));
            if ~isreal(ang(j)) % when normals are parallel
                ang(j) = 0;
            end
            fprintf(angleOut, '%f\n', ang(j)); 
        else
            ang(j) = 0;
            fprintf(angleOut, '%f\n', ang(j)); 
        end
        
        count = count + 1;
    end
        
end