clc;clear;close all;

tic
%Create a 3D domain
xmax = 300;
ymax = 300;
zmax = 200;
x = -2:1:xmax+3;  % Define x coordinates
y = -2:1:ymax+3;  % Define y coordinates
z = -2:1:zmax+3;  % Define z coordinates
dom_size = (length(x))*(length(y))*(length(y));
grid_points = zeros((dom_size),3);
RBLANK = ones(length(x),length(y),length(z));
% 0 solid
% 1 fluid
% 2 içteki yakın solid

c = 1;
for i=1:x(end)
    for j=1:y(end)
        for k=1:z(end)
            grid_points(c,:) = [x(i) y(j) z(k)];
            c = c +1;
        end
    end
end
toc
disp('Domain grid is created.')

%%
%contains every triangular surface mesh 3 points as cell array NX3
% tri_mesh_array = 31;
% tri_N = 100;
% tri_mesh_array_pts = (5)*rand(tri_N,9);
% norm_of = zeros(tri_N,1);
%
% tri_mesh_array(:,1) = mat2cell(tri_mesh_array_pts(:,1:3),[ones(1,tri_N)],[3]);
% tri_mesh_array(:,2) = mat2cell(tri_mesh_array_pts(:,4:6),[ones(1,tri_N)],[3]);
% tri_mesh_array(:,3) = mat2cell(tri_mesh_array_pts(:,7:9),[ones(1,tri_N)],[3]);
% tri_mesh_array(:,end+1) = num2cell(norm_of);

% load("tri_mesh_cell_naca6409.mat")
% load("vertices_naca6409.mat")
% load("tri_mesh_cell.mat") % loading tri_mesh_cell array contains each triangular surface mesh information
% %load("vertices.mat")
load('TRI_sphere_30dia.mat') % loading TRI structure that contains vertices and faces of solid
translational_vector = min(TRI.vertices(:,:));
TRI.vertices = TRI.vertices - translational_vector+[100 150 100];

% Tri_index_coord cell array generator
tri_index_coord=cell(1,4);
for i = 1:length(TRI.faces(:,1))
    indexes= TRI.faces(i,:);
    tri_index_coord(i,:)= {TRI.vertices(indexes(1),:), TRI.vertices(indexes(2),:), TRI.vertices(indexes(3),:), []};
end
tri_mesh_array = tri_index_coord;

% Pre allocation of necessary arrays
closest_solid_grid = zeros(length(tri_mesh_array),3);
closest_fluid_grid = zeros(length(tri_mesh_array),4);
NUP = zeros(18,length(tri_mesh_array));
NUP_TOT = zeros(length(tri_mesh_array),1);
DELTA = zeros(18,length(tri_mesh_array));
dir_vec = (1:18)';
search_R = 1;

% ref point is an arbitrary point inside the boundary solid shape
% this will be used for finding real triangular surface norm direction
ref_point = [mean(TRI.vertices(:,1)) mean(TRI.vertices(:,2)) mean(TRI.vertices(:,3))];

%finding each triangular mesh surface norrm direction

tic
c = 1;
cc = 1;
for i=1:length(tri_mesh_array)
    %tic

    % finding triangular mesh center
    tri_center_point = (tri_mesh_array{i,1}+ tri_mesh_array{i,2} + tri_mesh_array{i,3}) / 3;

    % finding triangular mesh normal vector 
    vec1 = tri_mesh_array{i,3} - tri_center_point;
    vec2 = tri_mesh_array{i,2} - tri_center_point;
    norm_vector = cross(vec1,vec2);
    tri_mesh_array{i,end} = norm_vector/norm(norm_vector);

    % correction of normal vector direction with using inner reference
    % point
    triangleToRefpoint_vec = ref_point - tri_center_point;
    triangleToRefpoint_vec = triangleToRefpoint_vec/norm(triangleToRefpoint_vec);
    % finding angle between candidate norm vector and reference vector
    angle = acos(dot(tri_mesh_array{i,end},triangleToRefpoint_vec));
    if angle < pi/2
        tri_mesh_array{i,end} = -tri_mesh_array{i,end};
    end

    % finding closest solid inner node with using triangular mesh normal
    % creating boundary cupe for every triangular mesh center point
    bound_max = floor(tri_center_point+search_R);
    bound_min = bound_max-search_R;
    boundcube = [
        bound_max;
        bound_max(1) bound_max(2) bound_max(3)-search_R  ;
        bound_max(1)-search_R bound_max(2) bound_max(3)  ;
        bound_max(1)-search_R bound_max(2) bound_max(3)-search_R;
        bound_min
        bound_min(1) bound_min(2) bound_min(3)+search_R  ;
        bound_min(1)+search_R bound_min(2) bound_min(3)  ;
        bound_min(1)+search_R bound_min(2) bound_min(3)+search_R
        ];

    % storing every point that is near to triangular
    %candidate_grid_point = grid_points(condition,:);
    candidate_grid_point= boundcube;

    % sepereta nearest inner and outer point relative triangular surface normal
    candidate_grid_point_rel_tri = candidate_grid_point - tri_center_point;
    % normalize the vectors
    candidate_grid_point_rel_tri_norm = candidate_grid_point_rel_tri ./ vecnorm(candidate_grid_point_rel_tri,2,2);

    % finding angle between candidate tri_to_grid vector and normal vector
    % dot product formula is used in below
    angle = acos(sum(tri_mesh_array{i,end} .* candidate_grid_point_rel_tri_norm,2));
    idx_inner = angle >= pi/2;

    %candidate_outer_grid_point_rel_tri = candidate_grid_point_rel_tri(idx_outer,:);
    candidate_inner_grid_point_rel_tri = candidate_grid_point_rel_tri(idx_inner,:);

    %candidate_outer_grid_point = candidate_grid_point(idx_outer,:);
    candidate_inner_grid_point = candidate_grid_point(idx_inner,:);

    % find nearest grid point to triangular
    [~ , idx_inner] = min(vecnorm(candidate_inner_grid_point_rel_tri,2,2));


    closest_solid_grid(i,:) = candidate_grid_point(idx_inner,:);

    % % Calculation of solid nearest nodes to fluid nodes
    direction_matrix = [ 1  0  0;   %right
                        -1  0  0;  %left
                         0 -1  0; %behind
                         0  1  0;   %front
                         0  0  1;   %top
                         0  0 -1;  %bottom
                         1 -1  0;  %right/behind
                        -1  1  0;  %left-front
                         1  1  0;  %right front
                        -1 -1  0; %left-behind
                         1  0  1;   %top right corner
                        -1  0 -1; %bottom left corner
                         1  0 -1 ; %bottom right corner
                        -1  0  1;  %top left corner
                         0 -1  1; %top behind corner
                         0  1 -1;  %bottom front corner
                         0 -1 -1; %bottom behind corner
                         0  1  1;   %top front corner
                                 ];
    direction_matrix_norm = direction_matrix ./ vecnorm(direction_matrix,2,2);
    candidate_fluid_nodes = closest_solid_grid(i,:) + direction_matrix;

    % calculation distance to fluid to boundary. used projection method
    support_vector = tri_center_point - closest_solid_grid(i,:);
    inner_to_bound_vector = sum(support_vector .* direction_matrix_norm,2).* direction_matrix_norm; % projection
    fluid_to_bound_dist = vecnorm(direction_matrix,2,2) - vecnorm(inner_to_bound_vector,2,2);
    angle = acos(sum(tri_mesh_array{i,end} .* direction_matrix_norm,2));
    idx = angle <= pi/2;
    fluid_to_bound_dist_org = fluid_to_bound_dist;
    fluid_to_bound_dist = fluid_to_bound_dist(idx);

    % find cloest fluid grid
    closest_fluid_grid(cc : cc + sum(idx) -1,:) = [candidate_fluid_nodes(idx,:) fluid_to_bound_dist];


    % NUP array generator
    dir_info = dir_vec .* idx;
    dir_info = dir_info(dir_info ~= 0);
    dir_col = zeros(18,1);
    dir_col(1:length(dir_info)) = dir_info;
    NUP(:,i) = dir_col;
    NUP_TOT(i) = sum(idx);
    % DELTA array generator
    DELTA(:,i) = fluid_to_bound_dist_org .* idx;
    [row_size, ~] = size(closest_solid_grid);
    NCURV = row_size;



    %idx = fluid_to_bound_dist < sqrt(2);

    %[~,idx] = min(fluid_to_bound_dist(idx));

    %closest_fluid_grid = [closest_fluid_grid ; [candidate_solid_nodes(idx,:) fluid_to_bound_dist]];
    %closest_fluid_grid(cc : cc + sum(idx) -1,:) = [candidate_fluid_nodes(idx,:) fluid_to_bound_dist];
    cc = cc+ sum(idx);
    %toc
end
toc
disp('Closest Inner Solid Nodes and Closest Fluid Nodes are found.')
IJK = closest_solid_grid';

%%
% Finding all solid nodes with using 'raycast' algorithm
% For performance efficiency we reduced grid search area for finding inner
% solid nodes to minimum cube that cover the whole solid shape
tic
cond1 = grid_points(:,1) > min(TRI.vertices(:,1)) & grid_points(:,1) < max(TRI.vertices(:,1));
cond2 = grid_points(:,2) > min(TRI.vertices(:,2)) & grid_points(:,2) < max(TRI.vertices(:,2));
cond3 = grid_points(:,3) > min(TRI.vertices(:,3)) & grid_points(:,3) < max(TRI.vertices(:,3));
condf = cond1 & cond2 & cond3;
grid_points_cube = grid_points(condf,:);
%max(TRI.vertices(:,1))
solid_grid_idx = in_polyhedron(TRI, grid_points_cube);
solid_grid_points = grid_points_cube(solid_grid_idx,:);

solid_grid_idx_3d = solid_grid_points + 3;
closest_solid_grid_idx_3d = closest_solid_grid + 3;

for i=1:length(solid_grid_idx_3d(:,1))
    RBLANK(solid_grid_idx_3d(i,1), solid_grid_idx_3d(i,2), solid_grid_idx_3d(i,3)) = 0;
%     RBLANK(closest_solid_grid_idx_3d(i,1), closest_solid_grid_idx_3d(i,1), closest_solid_grid_idx_3d(i,1)) = 2;
end

for i= 1:length(closest_solid_grid_idx_3d(:,1))
    RBLANK(closest_solid_grid_idx_3d(i,1), closest_solid_grid_idx_3d(i,2), closest_solid_grid_idx_3d(i,3))=2;
end
% boundary correction
for i=1:xmax
    for k=1:zmax
        RBLANK(i,ymax+1,k) = RBLANK(i,ymax,k);
        RBLANK(i,ymax+2,k) = RBLANK(i,ymax,k);
        RBLANK(i,ymax+3,k) = RBLANK(i,ymax,k);

        RBLANK(i,2,k) = RBLANK(i,1,k);
        RBLANK(i,1,k) = RBLANK(i,1,k);
    end
end


%RBLANK(ismember(grid_points,inner_grid_points,'rows')) = 0;

% RBLANK(ismember(grid_points,inner_grid_points,'rows')) = 0;
% RBLANK(~ismember(grid_points,inner_grid_points,'rows')) = 1;
% RBLANK(ismember(grid_points,closest_fluid_grid(:,1:3),'rows')) = 2;
toc
disp('All Domain Nodes is labeled as fluid/solid.')



%% Plot option
ds = 1000;
%grid_points_ds =  grid_points(1:ds:length(grid_points),1:ds:length(grid_points),1:ds:length(grid_points));
%plot3(grid_points_ds(:,1),grid_points_ds(:,2),grid_points_ds(:,3),'.',MarkerSize=1)
%hold on

% Plotting surface mesh vertices
plot3(TRI.vertices(:,1),TRI.vertices(:,2),TRI.vertices(:,3),'.g',MarkerSize=5)

% Plotting closest solid nodes
figure(1)
plot3(closest_solid_grid(:,1),closest_solid_grid(:,2),closest_solid_grid(:,3),'c.','MarkerSize',10)
daspect([1 1 1])
hold on 

% Plotting closest fluid nodes
% figure(2)
plot3(closest_fluid_grid(:,1),closest_fluid_grid(:,2),closest_fluid_grid(:,3),'m.','MarkerSize',10)
daspect([1 1 1])

% % Plotting for all inner solid nodes
% figure(1)
plot3(TRI.vertices(:,1), TRI.vertices(:,2),TRI.vertices(:,3),'.')
% hold on
% pt = inner_grid_points;
% daspect([1 1 1])
% plot3(pt(:,1), pt(:,2), pt(:,3),'.m')
xlim([0 300])
ylim([0 300])
zlim([0 200])

%% SAVE

save("IJK.mat","IJK")
save("RBLANK.mat","RBLANK")
save("DELTA.mat","DELTA")
save("D.mat","DELTA")
save("NCURV.mat","NCURV")
save("NUP.mat","NUP")
save("NUP_TOT.mat","NUP_TOT")


