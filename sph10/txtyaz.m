% % %% 10 Diameter Sphere In 50x50x50 Domain
% % 
% % load("RBLANK.mat") % Three dimensional array for each grid point. Contains zero for solid, one for fluid, two for solid grid points with fluid neighbours.
% % 
% % load("DELTA.mat")  % 18 rows for D3Q19 directions and columns represent solid grid points with fluid neighbours array contains fluid 
% %                    % grid point's distance to actual boundary. 
% % 
% % load("NUP.mat")    % Contains the neighbour directions for each boundary node
% % 
% % load("NUP_TOT.mat")% Neighbour count for each boundary node.
% % 
% % load("IJK.mat")    % Coordinates for each boundary node  
% % 
% % load("NCURV.mat")  % Boundary node count (integer)


%%Put space between indexes

% writematrix(RBLANK, 'RBLANK.txt', 'Delimiter', 'tab');
% writematrix(DELTA, 'D.txt', 'Delimiter', 'tab');
% writematrix(DELTA, 'DELTA.txt', 'Delimiter', 'tab');
% writematrix(NCURV, 'NCURV.txt', 'Delimiter', 'tab');
% writematrix(NUP, 'NUP.txt', 'Delimiter', 'tab');
% writematrix(NUP_TOT, 'NUP_TOT.txt', 'Delimiter', 'tab');
% writematrix(IJK, 'IJK.txt', 'Delimiter', 'tab');



RBLANK = RBLANK(:);
IJK = IJK(:);
DELTA = DELTA(:);
NUP = NUP(:);
NUP_TOT = NUP_TOT(:);



%%Put tab between indexes

writematrix(RBLANK, 'RBLANK.txt', 'Delimiter', ' ');
writematrix(DELTA, 'D.txt', 'Delimiter', ' ');
writematrix(DELTA, 'DELTA.txt', 'Delimiter', ' ');
% writematrix(NCURV, 'NCURV.txt', 'Delimiter', ' ');
writematrix(NUP, 'NUP.txt', 'Delimiter', ' ');
writematrix(NUP_TOT, 'NUP_TOT.txt', 'Delimiter', ' ');
writematrix(IJK, 'IJK.txt', 'Delimiter', ' ');
