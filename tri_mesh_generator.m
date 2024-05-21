

vertices = vertices;


min(vertices(:,:))

mean(vertices,1)
tri_index_coord=cell(1,4);

for i = 1:length(tri_index)
indexes=tri_index(i,:);
tri_index_coord(i,:)= {vertices(indexes(1),:), vertices(indexes(2),:), vertices(indexes(3),:), []};

end