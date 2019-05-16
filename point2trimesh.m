function [ distances, surface_points, faces2, vertices2, corresponding_vertices_ID, new_faces_ID ] = point2trimesh( varargin )
% -------------------------------------------------------
%
%    point2trimesh - Distance between a point and a triangulated surface in 3D 
%    
%    The shortest line connecting a point and a triangulation in 3D is 
%    computed. The nearest point on the surface as well as the distance
%    is returned. The distance is signed according to face normals to
%    identify on which side of the surface the query point resides. 
%    The implementation is optimized for speed, and depending on your
%    application you can use linear or parallel computation. 
%
%    Point insertion functionality 
%    (this feature is experimental and not optimized for speed):
%    If the function is called with more than two output arguments, 
%    the surface points are included into the given triangulation 
%    and Delaunay conditions are restored locally. If triangles with small
%    angles occur, additional vertices are inserted to eliminate them
%    if possible. 
%
%    Algorithm: From every query point, 
%       - the nearest vertex
%       - the nearest point on the edges and
%       - the nearest point on the triangle's surfaces
%    is calculated and the minimum distance out of these three is returned.   
%
%    Ver. 1.0.0
%
%    Created:         Daniel Frisch        (29.03.2015)
%    Last modified:   Daniel Frisch        (24.06.2015)
%
% ------------------------------------------------------
%
%  Inputs:
%      Pairs of parameter names and corresponding values. 
%      Structure arrays are expanded into separate inputs, 
%      where each field name corresponds to an input parameter name. 
%      Parameters:
%      - 'Faces'         (#faces    x 3) Triangulation connectivity list. Each row defines a face, elements are vertex IDs.  
%      - 'Vertices'      (#vertices x 3) Point matrix. Columns are x,y,z coordinates, row numbers are vertex IDs.  
%      - 'QueryPoints'   (#qPoints  x 3) Columns are x,y,z coordinates; each row defines a query point. Can be empty. 
%      - 'MaxDistance'   (1 x 1)         If the distance between a surface_point and its nearest vertex is within this range, 
%                                        no new vertex is inserted into the mesh. This helps avoiding
%                                        triangles with small angles. (default: 1/10 the smallest inradius)
%      - 'UseSubSurface' (1 x 1)         Logical. If true (default), the distance to edges and surfaces is only calculated
%                                        for faces that are connected to the vertex nearest to the query point.  
%                                        This speeds up the calculation but if the distance between two opposite parts 
%                                        of the surface is less than the spacing of the vertices, wrong results are produced.
%                                        In the vectorized algorithm, 'SubSurface' is always false.  
%      - 'Algorithm'     'linear' (default): query points are processed successively in a 'for' loop. 
%                            Use this if you have only few query points. 
%                        'parallel': query points are processed in parallel with 'parfor'. 
%                            If no parallel pool exists, Matlab creates one automatically (which takes half a minute).
%                            It shuts down after 30 min idle time by default, but you can change that in the "Parallel Preferences". 
%                            If Matlab doesn't correctly recognize the number of cores of your processor, 
%                            change the "Number of workers" in "Manage Cluster Profiles".   
%                         'vectorized': query points are processed altogether in a vectorized manner.
%                            Be careful, this needs much RAM if you have many query points. 
%                            This option is mostly included to show that vectorization does not always speed up the calculation:
%                            In the linear code, every query point can be assigned an individual cutout of the whole surface. 
%                         'linear_vectorized_subfunctions': query points and their individual cutout surfaces
%                            are processed successively in a for loop by functions that are capable of processing more than one point. 
%                            This option is included to show that the non-vectorized functions are faster 
%                            if only one point is processed at a time. 
%                         'parallel_vectorized_subfunctions': query points and their individual cutout surfaces
%                            are processed in parallel in a parfor loop by functions that are capable of processing more than one point. 
%                            Again, this option is included to show that the non-vectorized functions are faster 
%                            if only one point is processed at a time. 
%
%  Outputs:
%      - distances      (#qPoints   x 1)   Vector with the point-surface distances; sign depends on normal vectors. 
%      - surface_points (#qPoints   x 3)   Matrix with the corresponding nearest points on the surface. 
%      - faces2         (#faces2    x 3)   Connectivity matrix of the triangulation including the surface_points as dedicated vertices  
%      - vertices2      (#vertices2 x 3)   Point/Vertex matrix of the triangulation including the surface_points as dedicated vertices 
%      - corresponding_vertices_ID    (#qPoints x 1) Vector with the IDs of the vertices corresponding to the query points
%      - new_faces_ID                 Vector with the IDs of the new or modified faces (to give them a different color, for example)
%
%
% Usage example:
%      FV.faces    = [5 3 1; 3 2 1; 3 4 2; 4 6 2];
%      FV.vertices = [2.5 8.0 1; 6.5 8.0 2; 2.5 5.0 1; 6.5 5.0 0; 1.0 6.5 1; 8.0 6.5 1.5];
%      points      = [2 4 2; 4 6 2; 5 6 2];
%      [distances,surface_points] = point2trimesh(FV, 'QueryPoints', points); 
%      patch(FV,'FaceAlpha',.5); xlabel('x'); ylabel('y'); zlabel('z'); axis equal; hold on
%      plot3M = @(XYZ,varargin) plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),varargin{:});
%      plot3M(points,'*r')
%      plot3M(surface_points,'*k')
%      plot3M(reshape([shiftdim(points,-1);shiftdim(surface_points,-1);shiftdim(points,-1)*NaN],[],3),'k')
%
%
% Dependencies
% [flist,plist] = matlab.codetools.requiredFilesAndProducts('point2trimesh.m'); [flist'; {plist.Name}']
%
% (no dependencies) 
%
% Parse Inputs 
% valdiation functions
faceChk   = @(x) validateattributes(x,{'double','int32' },{'real','finite','nonnan','positive'   ,'integer','size',[NaN 3]});
vertChk   = @(x) validateattributes(x,{'double','single'},{'real','finite','nonnan'                        ,'size',[NaN 3]});
pointChk  = @(x) validateattributes(x,{'double'         },{'real','finite','nonnan'                                       });
distChk   = @(x) validateattributes(x,{'double'         },{'real','finite','nonnan','nonnegative'          ,'scalar'      });
logicChk  = @(x) validateattributes(x,{'logical'        },{'scalar'});
charChk   = @(x) validateattributes(x,{'char'           },{'nonempty','vector'});
parser = inputParser;
parser.FunctionName = mfilename;
parser.addParameter('Faces'         ,[]      ,faceChk);
parser.addParameter('Vertices'      ,[]      ,vertChk);
parser.addParameter('QueryPoints'   ,[]      ,pointChk);
parser.addParameter('MaxDistance'   ,[]      ,distChk);
parser.addParameter('UseSubSurface' ,true    ,logicChk);
parser.addParameter('Algorithm'     ,'linear',charChk);
parser.parse(varargin{:});
faces    = double(parser.Results.Faces);
vertices = double(parser.Results.Vertices);
assert(~isempty(faces) && ~isempty(vertices), 'Invalid argument: ''Faces'' and ''Vertices'' mustn''t be empty.')
assert(max(faces(:))<=size(vertices,1), 'The value of ''Faces'' is invalid: the maximum vertex ID is bigger than the number of vertices in ''Vertices''')
qPoints = parser.Results.QueryPoints;
useSubSurface = parser.Results.UseSubSurface;

if nargout>2, insertPoints = true;
else
    insertPoints = false; 
end
% Calculate normals
r1 = vertices(faces(:,1),:);  % (#faces x 3) % 1st vertex of every face
r2 = vertices(faces(:,2),:);  % (#faces x 3) % 2nd vertex of every face
r3 = vertices(faces(:,3),:);  % (#faces x 3) % 3rd vertex of every face
normals = cross((r2-r1),(r3-r1),2); % (#faces x 3) normal vector of every face
normals = bsxfun(@rdivide,normals,sqrt(sum(normals.^2,2))); % (#faces x 3) normalized normal vector of every face
if isempty(qPoints)
    distances = [];
    surface_points = [];
    faces2 = faces;
    vertices2 = vertices;
    corresponding_vertices_ID = [];
    new_faces_ID = [];
    return
end
% Distance Calculation
nQPoints = size(qPoints,1);
D = NaN(nQPoints,1);
P = NaN(nQPoints,3);
if insertPoints
    max_distance = parser.Results.MaxDistance;
    if isempty(max_distance)
        tri = triangulation(faces,vertices);
        [~,r] = tri.incenter;
        max_distance = min(r)/10;
    end
    max_distance = max(max_distance,100*eps);
    is_new_face = zeros(size(faces,1),1);
    is_new_vertex = zeros(size(vertices,1),1);
    corresponding_vertices_ID = NaN(size(qPoints,1),1);
end
switch parser.Results.Algorithm
    case {'linear','normal'}
        for r = 1:nQPoints
            % Determine the surface points
            point = qPoints(r,:); % (1 x 3) query point
            [d,p,f] = processPoint(faces,vertices,point,normals, @distance_to_vertices,@distance_to_edges,@distance_to_surfaces, useSubSurface);
            D(r) = d;
            P(r,:) = p;
            if insertPoints
                % Include the surface points as new vertices into the mesh and restore Delaunay conditions 
                [ tri, is_new_face, is_new_vertex, new_vertex_ID ] = insert_vertex( triangulation(faces,vertices), p, f, max_distance, is_new_face, is_new_vertex );
                faces    = tri.ConnectivityList;
                vertices = tri.Points;  
                normals  = tri.faceNormal;
                corresponding_vertices_ID(r) = new_vertex_ID;
            end
        end
        
    case 'parallel'
        assert(~insertPoints,'''Algorithm'', ''%s'' doesn''t support including the surface points into the geometry. Call point2trimesh with fewer output arguments or use ''linear'' algorithm.',parser.Results.Algorithm)
        parfor r = 1:nQPoints
            point = qPoints(r,:); % (1 x 3) query point
            [d,p] = processPoint(faces,vertices,point,normals, @distance_to_vertices,@distance_to_edges,@distance_to_surfaces, useSubSurface);
            D(r) = d;
            P(r,:) = p;
        end
        
    case 'vectorized'
        assert(~insertPoints,'''Algorithm'', ''%s'' doesn''t support including the surface points into the geometry. Call point2trimesh with fewer output arguments or use ''linear'' algorithm.',parser.Results.Algorithm)
        if useSubSurface && ~ismember('UseSubSurface',parser.UsingDefaults)
            warning('You specified ''UseSubSurface'',true, but ''Algorithm'',''vectorized'' always searches on the complete surface') 
        end
        [D1,P1] = distance_to_vertices_vectorized(faces,vertices,qPoints,normals); % (#qPoints x 1), (#qPoints x 3), (#qPoints x 1) 
        [D2,P2] = distance_to_edges_vectorized   (faces,vertices,qPoints,normals);
        [D3,P3] = distance_to_surfaces_vectorized(faces,vertices,qPoints,normals);
        % find minimum distance type
        D = [D1,D2,D3];      % (#qPoints x 3)
        P = cat(3,P1,P2,P3); % (#qPoints x xyz x 3)
        [~,I] = min(abs(D),[],2);
        D = D(sub2ind(size(D),(1:length(I))',I));
        % extract nearest point on surface
        P = permute(P,[2,3,1]); % (xyz x 3 x #qPoints)
        sz = [size(P) 1 1];
        P = P(:,sub2ind(sz(2:3),I,(1:length(I))')); % (xyz x #qPoints)
        P = P'; % (#qPoints x xyz)
        
    case 'linear_vectorized_subfunctions'
        for r = 1:nQPoints
            % Determine the surface points
            point = qPoints(r,:); % (1 x 3) query point
            [d,p,f] = processPoint(faces,vertices,point,normals, @distance_to_vertices_vectorized,@distance_to_edges_vectorized,@distance_to_surfaces_vectorized, useSubSurface);
            D(r) = d;
            P(r,:) = p;
            if insertPoints
                % Include the surface points as new vertices into the mesh and restore Delaunay conditions 
                [ tri, is_new_face, is_new_vertex, new_vertex_ID ] = insert_vertex( triangulation(faces,vertices), p, f, max_distance, is_new_face, is_new_vertex );
                faces    = tri.ConnectivityList;
                vertices = tri.Points;  
                normals  = tri.faceNormal;
                corresponding_vertices_ID(r) = new_vertex_ID;
            end
        end
        
    case 'parallel_vectorized_subfunctions'
        assert(~insertPoints,'''Algorithm'', ''%s'' doesn''t support including the surface points into the geometry. Call point2trimesh with fewer output arguments or use ''linear'' algorithm.',parser.Results.Algorithm)
        parfor r = 1:nQPoints
            point = qPoints(r,:); % (1 x 3) query point
            [d,p] = processPoint(faces,vertices,point,normals, @distance_to_vertices_vectorized,@distance_to_edges_vectorized,@distance_to_surfaces_vectorized, useSubSurface);
            D(r) = d;
            P(r,:) = p;
        end
        
    otherwise
        error('The value of ''Algorithm'' is invalid.')
end
if insertPoints
    % Despite the Delaunay condition is fulfilled, the
    % triangulation might still contain triangles with small angles.
    % To prevent this, insert vertices at the circumcenters  of these triangles.
    tri = triangulation(faces,vertices);
    while true
        [minAng,worstFace] = minimumAngle(tri,(1:size(tri.ConnectivityList,1))');
        insertPt = tri.circumcenter(worstFace);
        [~,insertPt,f] = processPoint(tri.ConnectivityList,tri.Points,insertPt,tri.faceNormal, @distance_to_vertices,@distance_to_edges,@distance_to_surfaces, useSubSurface);        
        [~,r] = tri.incenter(f);
        [ tri2, is_new_face2, is_new_vertex2 ] = insert_vertex( tri, insertPt, f, r/2, is_new_face, is_new_vertex );
        [minAng2,~] = minimumAngle(tri2,(1:size(tri2.ConnectivityList,1))');
        if minAng2 > minAng+eps
            tri = tri2;
            is_new_face   = is_new_face2;
            is_new_vertex = is_new_vertex2;
        else
            break;
        end
    end
    faces2    = tri.ConnectivityList;
    vertices2 = tri.Points;
    new_faces_ID = find(is_new_face);
    
end
% return output arguments
distances      = D;  % (#qPoints x 1)
surface_points = P;  % (#qPoints x 3)
end
