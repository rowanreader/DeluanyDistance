function run
file = 'sphere_points_leftBreast_S1.txt';
main(file); % runs vectors and initialize dirNum
results('output.txt'); % computes surface points
%angles();
circleFit(); % uses surface points to fit a sphere and get new dots
results('output2.txt'); % gets new surface points
angles(); % computes angles
