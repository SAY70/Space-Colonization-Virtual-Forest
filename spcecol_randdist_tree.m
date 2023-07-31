clearvars; close all; clc;

num_trees = 100;  % number of trees to generate
b_m      = 0.7;  % branch magnitude (for reach node segment)
min_d    =  1;   % minimum distance to branch
max_d    =  2;   % maximum distance to branch
max_its  = 20;
min_h = 13; % minimum height of the tree
max_h = 16; % maximum height of the tree
min_buds = 80; % minimum numbers of tree buds
max_buds = 120; % maximum numbers of tree buds
cone_r = 1.7; % radius of the cone (bottom r)
% top_r = 0.5; % radius of the cone at the top
mid_r = 2.5; % radius at the middle of the cone
bottom_r = 1.5; % Note: This is top_r /need to replace
bud_start_pos=4;

% Create a figure and plot the trees
figure('color', 'white');
hold on;
%% 
%  x_offsets = [-5, 6, 8, 0, 5, -3];
%  y_offsets = [-9, -9, 3, -2, -3, 4];
%z_offsets = 0; 

rng(123774); % set seed value
%%
% Defining number of cells in the x and y directions
num_cells_x = 10;
num_cells_y = 10;
max_dist=0.5;
min_dist=0.5;

% Calculating size of each cell
cell_size_x = round(100 / num_cells_x);
cell_size_y = round(100 / num_cells_y);

% Generating random offsets for each cell
cell_offsets_x = rand(num_trees, 1) * (cell_size_x - (max_dist + min_dist)) + min_dist;
cell_offsets_y = rand(num_trees, 1) * (cell_size_y - (max_dist + min_dist)) + min_dist;

%%
% Initialize arrays to store all head and tail positions
all_head_positions = [];
all_tail_positions = [];

%%
trunk_head_pos = [];
trunk_tail_pos = [];
primary_head_pos = [];
primary_tail_pos = [];
secondary_head_pos = [];
secondary_tail_pos = [];
tertiary_head_pos = [];
tertiary_tail_pos = [];
no_use_head_pos = [];
no_use_tail_pos = [];

%%
%elevation angles
trunk_el_angles=[];
primary_el_angles=[];
secondary_el_angles=[];
tertiary_el_angles=[];
no_use_el_angles=[];

%%
trunk_centers = [];
primary_centers = [];
secondary_centers = [];
tertiary_centers = [];
no_use_centers = [];

%%
all_radius = [];
trunk_radius=[];
primary_radius = [];
secondary_radius = [];
tertiary_radius = [];
no_use_radius = [];

all_branch_lengths=[];
trunk_lengths=[];
primary_lengths = [];
secondary_lengths = [];
tertiary_lengths = [];
no_use_lengths = [];

cone_h=[];
for i = 1:num_trees
    cone_h = round((max_h - min_h) * rand + min_h);
    num_buds = round(max_buds - (cone_h/max_h) * (max_buds - min_buds));
    theta = 2*pi*rand(num_buds,1);
%     r = cone_r*sqrt(rand(num_buds,1));
%     h = cone_h*rand(num_buds,1);
%     taper = (h/cone_h) .* (top_r-cone_r)/cone_h; % taper the radius towards the top
%     tapered_r = cone_r + taper.*h + (1-taper).*r;

    r = cone_r*sqrt(rand(num_buds,1));
    h = cone_h*rand(num_buds,1);
    taper = ((h < cone_h/2) .* (mid_r-cone_r)/cone_h + (h >= cone_h/2) .* (bottom_r-mid_r)/cone_h); % taper the radius towards the top
    tapered_r = cone_r + taper.*h + (1-taper).*r;

    b_x = tapered_r.*cos(theta);
    b_y = tapered_r.*sin(theta);
    b_z = h+bud_start_pos;

    % Create SC_tree object and set the inputs
    sc_tree(i) = SC_tree_obj();
    sc_tree(i).set_inputs([b_x'; b_y'; b_z'], b_m, min_d, max_d, max_its);
    sc_tree(i).make_root;
    sc_tree(i).make_branch;
    sc_tree(i).consolidate_vectors;
    sc_tree(i).identify_branches;
    sc_tree(i).shinozaki_pipe(4,0.9);
    sc_tree(i).assign_diel();
    sc_tree(i).get_pos();
    vol = sc_tree(i).get_vol();
    bh = sc_tree(i).tree.heads;
    bt = sc_tree(i).tree.tails;
    as = sc_tree(i).semi_major;
    plt3k = @(a,w) plot3(a(1,:), a(2,:), a(3,:),'color', 'k', 'linewidth', w);
    
    % Randomly position the tree in the box
%     x_offset = 30*rand(1)-10;
%     y_offset = 30*rand(1)-10;
%     z_offset = 0;

%     x_offset = linspace(1.5, 4, num_trees);
%     y_offset = linspace(1.5, 4, num_trees);
%     z_offset = 0;
% 
%       x_offset = x_offsets(i);
%       y_offset = y_offsets(i);
%       z_offset = 0; 
% Calculate position of tree in cell
    cell_x = mod(i-1, num_cells_x) + 1;
    cell_y = ceil(i/num_cells_x);
    x_offset = (cell_x - 1) * cell_size_x + cell_offsets_x(i);
    y_offset = (cell_y - 1) * cell_size_y + cell_offsets_y(i);
    z_offset = 0;
   
    for j = 1:size(bh,2)
        cr = as(j);
        ch = bh(:,j)+ [x_offset;y_offset;z_offset];
        ct = bt(:,j)+ [x_offset;y_offset;z_offset];
        all_head_positions = [all_head_positions ch];
        all_tail_positions = [all_tail_positions ct];
        all_radius = [all_radius 3*cr];
        plt3k([ch ct], cr);
        hold on; grid on;
        branch_lengths = norm(ch - ct);
        center_pos = ct + (ch - ct) / 2;
        z_axis = [0 0 1];
        all_branch_lengths=[all_branch_lengths branch_lengths];
        if cr > 3.3 % trunk branch
        trunk_radius = [trunk_radius 3*cr]; %for pine tree scaling factor 3
        trunk_lengths = [trunk_lengths branch_lengths];  
        trunk_centers = [trunk_centers, center_pos];
        trunk_heads = ch; trunk_tails = ct;
        trunk_head_pos = [trunk_head_pos trunk_heads];
        trunk_tail_pos = [trunk_tail_pos trunk_tails];
        trunk_dir = trunk_heads - trunk_tails;
        cos_elevation = dot(z_axis, trunk_dir) / norm(trunk_dir);
        trunk_elevation_deg = (acos(cos_elevation))* 180 / pi;
        trunk_el_angles = [trunk_el_angles trunk_elevation_deg];
        else
            if (cr >= 1.9) && (cr <= 3.3)% primary branch
            primary_radius = [primary_radius 3*cr];
            primary_lengths = [primary_lengths branch_lengths];
            primary_centers = [primary_centers, center_pos];
            primary_heads = ch;
            primary_tails = ct;
            primary_head_pos = [primary_head_pos primary_heads];
            primary_tail_pos = [primary_tail_pos  primary_tails];
            primary_dir = primary_heads - primary_tails;
            primary_cos_elevation = dot(z_axis, primary_dir) / norm(primary_dir);
            primary_elevation_rad = acos(primary_cos_elevation);
            primary_elevation_deg = primary_elevation_rad * 180 / pi;
            primary_el_angles = [primary_el_angles primary_elevation_deg];
            else
                if (cr >= 1.4) && (cr <= 1.9)% secondary branch
                secondary_radius = [secondary_radius 3*cr];
                secondary_lengths = [secondary_lengths branch_lengths];
                secondary_centers = [secondary_centers, center_pos];
                secondary_heads = ch;
                secondary_tails = ct;
                secondary_head_pos = [secondary_head_pos secondary_heads];
                secondary_tail_pos = [secondary_tail_pos  secondary_tails];
                secondary_dir = secondary_heads - secondary_tails;
                secondary_cos_elevation = dot(z_axis, secondary_dir) / norm(secondary_dir);
                secondary_elevation_rad = acos(secondary_cos_elevation);
                secondary_elevation_deg = secondary_elevation_rad * 180 / pi;
                secondary_el_angles = [secondary_el_angles secondary_elevation_deg];
                else
                    if (cr >= 1.1) && (cr <= 1.4)% tertiary branch
                tertiary_radius = [tertiary_radius 3*cr];
                tertiary_lengths = [tertiary_lengths branch_lengths];
                tertiary_centers = [tertiary_centers, center_pos];
                tertiary_heads = ch;
                tertiary_tails = ct;
                tertiary_head_pos = [tertiary_head_pos tertiary_heads];
                tertiary_tail_pos = [tertiary_tail_pos  tertiary_tails];
                tertiary_dir = tertiary_heads - tertiary_tails;
                tertiary_cos_elevation = dot(z_axis, tertiary_dir) / norm(tertiary_dir);
                tertiary_elevation_rad = acos(tertiary_cos_elevation);
                tertiary_elevation_deg = tertiary_elevation_rad * 180 / pi;
                tertiary_el_angles = [tertiary_el_angles tertiary_elevation_deg];
                else
                no_use_radius = [no_use_radius 3*cr];
                no_use_lengths = [no_use_lengths branch_lengths];
                no_use_centers = [no_use_centers, center_pos];
                no_use_heads = ch;
                no_use_tails = ct;
                no_use_head_pos = [no_use_head_pos no_use_heads];
                no_use_tail_pos = [no_use_tail_pos  no_use_tails];
                no_use_dir = no_use_heads - no_use_tails;
                no_use_cos_elevation = dot(z_axis, no_use_dir) / norm(no_use_dir);
                no_use_elevation_rad = acos(no_use_cos_elevation);
                no_use_elevation_deg = no_use_elevation_rad * 180 / pi;
                no_use_el_angles = [no_use_el_angles no_use_elevation_deg];
                end
            end
        end
    end
end

    view(3)
    zlim([0, cone_h+bud_start_pos]);
    xlabel('in Meter'); ylabel('in Meter'); zlabel('in Meter');
    axis equal;
    
end

%%
TPS_center_points = [trunk_centers, primary_centers, secondary_centers, tertiary_centers, no_use_centers];
%scatter3(TPS_center_points(1,:), TPS_center_points(2,:), TPS_center_points(3,:),'b.')
hold on
%%
receiver_positions = [];
receiver_positions_manual = [[50	85	0.6]; 
    [61.3644814221639	83.1036034595222	0.6]; 
    [71.4974449441384	77.6199178288738	0.6]; 
    [79.3008267391885	69.1431855342849	0.6];  
    [83.9290093078766	58.5919920499280	0.6]; 
    [84.8804572552334	47.1097229084684	0.6]; 
    [82.0520664329270	35.9406601371461	0.6]; 
    [75.7503368735596	26.2951449930991	0.6]; 
    [66.6581587562976	19.2184187077729	0.6]; 
    [55.7608106598257	15.4773543809047	0.6]; 
    [44.2391893401743	15.4773543809047	0.6]; 
    [33.3418412437024	19.2184187077729	0.6]; 
    [24.2496631264404	26.2951449930991	0.6]; 
    [17.9479335670730	35.9406601371461	0.6]; 
    [15.1195427447666	47.1097229084684	0.6]; 
    [16.0709906921234	58.5919920499280	0.6];
    [20.6991732608115	69.1431855342849	0.6];
    [28.5025550558616	77.6199178288738	0.6];
    [38.6355185778361	83.1036034595222	0.6];
    [50.0000000000000	85	0];
];
%n_scans = 1;

% for i = 1:n_scans
%     for j = 1:size(receiver_positions_manual, 1)
%         receiver_position = receiver_positions_manual(j,:);
%         receiver_positions = [receiver_positions; receiver_position];
%         receiver_position = receiver_position + [12, 0, 0];
%     end
% end



%%
% receiver_positions= [];
% receiver_position = [-5, 20, 0.6];
% n_scans = 3;

%for i = 1:n_scans
%     receiver_position = receiver_position + [12, 0, 0]; % Example: move 0.1 units to the right in the x direction
%     receiver_positions = [receiver_positions; receiver_position];

%     for j = 1:size(receiver_positions, 1)
%         receiver_position = receiver_positions(j,:);
    for j = 1:size(receiver_positions_manual, 1)
        receiver_position = receiver_positions_manual(j,:);
        %scatter3(receiver_position(1,1), receiver_position(1,2), receiver_position(1,3),'filled','ro');
        text(receiver_position(1), receiver_position(2), receiver_position(3)+0.2, num2str(j), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'r', 'FontSize', 18, 'FontWeight', 'bold');
        if j > 1
            previous_position = receiver_positions_manual(j-1,:);
            line([previous_position(1), receiver_position(1)], [previous_position(2), receiver_position(2)], [previous_position(3), receiver_position(3)], 'Color', 'r');
        end
        view_angle = 20;
    
        % Generate L-fractal forest particle center position
        particle_positions = TPS_center_points';
        
        % Calculate particle vectors relative to receiver location
        particle_vectors = particle_positions - receiver_position;
        
        % Calculate angles between particle vectors and receiver view direction
        receiver_direction = [0, 0, 1]; % upward direction
        cos_angles = dot(particle_vectors, repmat(receiver_direction, size(particle_vectors, 1), 1), 2) ./ vecnorm(particle_vectors, 2, 2);
        angles = acosd(cos_angles);
        
        
        % Define cylindrical surface parameters
        cyl_height = max(particle_positions(:,3)); % height of the cylinder
        cyl_radius = cyl_height * tand(view_angle); % radius of the cylinder
        
        
        % Create the cylindrical surface
        n_points = 20;
        theta = linspace(0, 2*pi, n_points);
        [X,Y,Z] = cylinder(cyl_radius, n_points);
        Z = Z * cyl_height + receiver_position(3);
        
        % Plot the particles in the cylindrical surface
        in_cyl = find((particle_positions(:,1)-receiver_position(1)).^2 + (particle_positions(:,2)-receiver_position(2)).^2 <= cyl_radius^2 ...
            & particle_positions(:,3) >= receiver_position(3) & particle_positions(:,3) <= receiver_position(3)+cyl_height);
        
        scatter3(particle_positions(in_cyl,1), particle_positions(in_cyl,2), particle_positions(in_cyl,3), 'filled', 'g');
        hold on;
        cyl_surface = surf(X+receiver_position(1), Y+receiver_position(2), Z, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        title('Particles in the receiver view and within the cylindrical surface');
        xlabel('In Meters');
        ylabel('In Meters');
        zlabel('In Meters');
        axis equal;
        grid on;
    
        trunk_in_view = find(ismember(trunk_centers', ...
            TPS_center_points(:,in_cyl)', 'rows'));
        primary_in_view = find(ismember(primary_centers', ...
            TPS_center_points(:,in_cyl)', 'rows'));
        secondary_in_view = find(ismember(secondary_centers', ...
            TPS_center_points(:,in_cyl)', 'rows'));
        tertiary_in_view = find(ismember(tertiary_centers', ...
            TPS_center_points(:,in_cyl)', 'rows'));
        no_use_in_view = find(ismember(no_use_centers', ...
            TPS_center_points(:,in_cyl)', 'rows'));
        
        % Find the positions of the primary and secondary vectors
        trunk_in_view_pos = particle_positions(trunk_in_view,:);
        primary_in_view_pos = particle_positions(primary_in_view, :);
        secondary_in_view_pos = particle_positions(secondary_in_view, :);
        tertiary_in_view_pos = particle_positions(tertiary_in_view, :);
        no_use_in_view_pos = particle_positions(no_use_in_view, :);
        
        % Calculate the volume of the cylinder
        cyl_volume = pi * cyl_radius^2 * cyl_height;
        
        % Calculate the number density of each particle type within the cylindrical surface
        % Trunk att
        trunk_num_density = (length(trunk_in_view) / (pi * cyl_radius^2));
        T1.DNSTY(j) = trunk_num_density;
        trunk_mean_length = mean(trunk_lengths(trunk_in_view));
        T1.LEN(j) = trunk_mean_length;
        trunk_mean_radius = 10^-2*mean(trunk_radius(trunk_in_view));
        T1.RAD(j) = trunk_mean_radius;
        min_el_trunk= min(trunk_el_angles(trunk_in_view));
        T1.MIN_EL(j) = min_el_trunk;
        max_el_trunk= max(trunk_el_angles(trunk_in_view));
        T1.MAX_EL(j) = max_el_trunk;
        
        % Primary att
        primary_num_density = (length(primary_in_view) / cyl_volume);
        B1.DNSTY(j) = primary_num_density;
        primary_mean_length = mean(primary_lengths(primary_in_view));
        B1.LEN(j) = primary_mean_length;
        primary_mean_radius = 10^-2*mean(primary_radius(primary_in_view));
        B1.RAD(j) = primary_mean_radius;
        min_el_primary= min(primary_el_angles(primary_in_view));
        B1.MIN_EL(j) = min_el_primary;
        max_el_primary= max(primary_el_angles(primary_in_view));
        B1.MAX_EL(j) = max_el_primary;
        
        %SEC att
        secondary_num_density = (length(secondary_in_view) / cyl_volume);
        B2.DNSTY(j) = secondary_num_density;
        secondary_mean_length = mean(secondary_lengths(secondary_in_view));
        B2.LEN(j) = secondary_mean_length;
        secondary_mean_radius = 10^-2*mean(secondary_radius(secondary_in_view));
        B2.RAD(j) = secondary_mean_radius;
        min_el_secondary= min(secondary_el_angles(secondary_in_view));
        B2.MIN_EL(j) = min_el_secondary;
        max_el_secondary= max(secondary_el_angles(secondary_in_view));
        B2.MAX_EL(j) = max_el_secondary;
        
        %Ter Att
        tertiary_num_density = (length(tertiary_in_view ) / cyl_volume);
        B3.DNSTY(j) = tertiary_num_density;
        tertiary_mean_length = mean(tertiary_lengths(tertiary_in_view));
        B3.LEN(j) = tertiary_mean_length;
        tertiary_mean_radius = mean(tertiary_radius(tertiary_in_view));
        B3.RAD(j) = 10^-2*tertiary_mean_radius;
        min_el_tertiary= min(tertiary_el_angles(tertiary_in_view));
        B3.MIN_EL(j) = min_el_tertiary;
        max_el_tertiary= max(tertiary_el_angles(tertiary_in_view));
        B3.MAX_EL(j) = max_el_tertiary;
        
        %No use att
        no_use_num_density = (length(no_use_in_view) / cyl_volume);
        B4.DNSTY(j) = no_use_num_density;
        no_use_mean_length = mean(no_use_lengths(no_use_in_view));
        B4.LEN(j) = tertiary_mean_length;
        no_use_mean_radius = 10^-2*mean(no_use_radius(no_use_in_view));
        B4.RAD(j) = no_use_mean_radius;
        min_el_no_use= min(no_use_el_angles(no_use_in_view));
        B4.MIN_EL(j) = min_el_no_use;
        max_el_no_use= max(no_use_el_angles(no_use_in_view));
        B4.MAX_EL(j) = max_el_no_use;
    end
%end
%%
% % Create a table with the results
% table = array2table([T1.DNSTY',T1.LEN',T1.RAD',T1.MIN_EL',T1.MAX_EL', ...
%     B1.DNSTY',B1.LEN',B1.RAD',B1.MIN_EL',B1.MAX_EL', ...
%     B2.DNSTY',B2.LEN',B2.RAD',B2.MIN_EL',B2.MAX_EL', ...
%     B3.DNSTY',B3.LEN',B3.RAD',B3.MIN_EL',B3.MAX_EL', ...
%     B4.DNSTY',B4.LEN',B4.RAD',B4.MIN_EL',B4.MAX_EL'], ...
%     'VariableNames', {'T1.DNSTY','T1.LEN','T1.RAD','T1.MIN_EL','T1.MAX_EL', ...
%     'B1.DNSTY','B1.LEN','B1.RAD','B1.MIN_EL','B1.MAX_EL', ...
%     'B2.DNSTY','B2.LEN','B2.RAD','B2.MIN_EL','B2.MAX_EL', ...
%     'B3.DNSTY','B3.LEN','B3.RAD','B3.MIN_EL','B3.MAX_EL', ...
%     'B4.DNSTY','B4.LEN','B4.RAD','B4.MIN_EL','B4.MAX_EL'}, 'RowNames', {'Receiver1', 'Receiver2', 'Receiver3','Receiver4','Receiver5','Receiver6','Receiver7','Receiver8','Receiver9','Receiver10','Receiver11','Receiver12','Receiver13','Receiver14','Receiver15','Receiver16','Receiver17','Receiver18','Receiver19','Receiver20'});
% 
% % Specify the file name and path to save the Excel file
% filename = 'forest_attributes_space.xlsx';
% path = 'C:\Users\say70\OneDrive - Mississippi State University\Desktop\PIERS\scobi_og\SC_tree_stat\';
% full_path = fullfile(path, filename);
% 
% % Write the table to an Excel file
% writetable(table, full_path, 'Sheet', 'Sheet1','WriteRowNames', 1);
% 
% 
% 
% % histogram
% % Initialize variables to store the radius and length data
% all_radius = [trunk_radius, primary_radius, secondary_radius, tertiary_radius, no_use_radius];
% all_branch_lengths = [trunk_lengths, primary_lengths, secondary_lengths, tertiary_lengths, no_use_lengths];
% 
% % Set colors and labels for the different histograms
% color_trunk='y';
% color_primary = 'r';
% color_secondary = 'g';
% color_tertiary = 'b';
% color_no_use = 'm';
% %legend_labels = {'Trunk','Primary Branches', 'Secondary Branches', 'Tertiary Branches', 'Qauternary Branches'};
% 
% % Plot the histograms for radius and length
% figure;
% subplot(2,1,1);
% %histogram(trunk_radius, 'BinWidth', 0.1, 'FaceColor', color_trunk);
% legend_labels = {'Qauternary Branches'};
% % hold on;
% %histogram(primary_radius, 'BinWidth', 0.1, 'FaceColor', color_primary);
% %histogram(secondary_radius, 'BinWidth', 0.1, 'FaceColor', color_secondary);
%  %histogram(tertiary_radius, 'BinWidth', 0.1, 'FaceColor', color_tertiary);
% histogram(no_use_radius, 'BinWidth', 0.1, 'FaceColor', color_no_use);
% title('Radius Histogram');
% xlabel('Radius (cm)');
% ylabel('Frequency');
% legend(legend_labels);
% grid on;
% 
% subplot(2,1,2);
% figure;
% %histogram(trunk_lengths, 'BinWidth', 0.1, 'FaceColor', color_trunk);
% hold on;
% legend_labels = {'Quaternary Branches'};
% %histogram(primary_lengths, 'BinWidth', 0.1, 'FaceColor', color_primary);
% %histogram(secondary_lengths, 'BinWidth', 0.1, 'FaceColor', color_secondary);
% histogram(tertiary_lengths, 'BinWidth', 0.1, 'FaceColor', color_tertiary);
% %histogram(no_use_lengths, 'BinWidth', 0.1, 'FaceColor', color_no_use);
% title('Branch Length Histogram');
% xlabel('Branch Length (m)');
% ylabel('Frequency');
% legend(legend_labels);
% grid on;
