function [node_data] = assign_pop_2_nodes(node_data, boundary, pop_Data, pop_R)
%ASSIGN_POP_2_NODES assign the population to the road nodes, based on 
% a raster population data
%   Input: node_data: the node data, without the population data
%          boundary: the boundary of the city
%          pop_Data: the population data
%          pop_R: the reference of the population data
%   Output: node_data: the node data, with the population data

% 20230305: add the judgment whether the cell is within the boundary;

% Set negative population values to zero
pop_Data(pop_Data < 0) = 0;

% Get the latitude and longitude steps from the raster reference
pop_lat_step = pop_R.CellExtentInLatitude;
pop_log_step = pop_R.CellExtentInLongitude;

% Determine the bounding box of the boundary
min_lat = min(boundary.Y);
max_lat = max(boundary.Y);
min_log = min(boundary.X);
max_log = max(boundary.X);

% Find the index range in pop_Data that corresponds to the bounding box
min_yid = max(1, floor((pop_R.LatitudeLimits(2) - max_lat) / pop_lat_step) + 1);
max_yid = min(size(pop_Data, 1), floor((pop_R.LatitudeLimits(2) - min_lat) / pop_lat_step) + 1);
min_xid = max(1, floor((min_log - pop_R.LongitudeLimits(1)) / pop_log_step) + 1);
max_xid = min(size(pop_Data, 2), floor((max_log - pop_R.LongitudeLimits(1)) / pop_log_step) + 1);

% Initialize the fourth column of node_data to zero (for population data)
node_data(:, 4) = 0;
sample_step=10;

init_node_id = node_data(:, 1);
node_data(:,1) = 1:length(node_data(:,1));

% Iterate over the filtered subset of pop_Data
for yid = min_yid:max_yid
    % Update progress display
    progress_percent = ((yid-min_yid+1) / (max_yid-min_yid+1)) * 100;
    % Construct the progress message
    % fprintf('\rProcessed percentage %.2f%%', progress_percent);

    for xid = min_xid:max_xid

        
        if ~isnan(pop_Data(yid, xid))
            [sample_lat, sample_log] = intrinsicToGeographic(pop_R, xid, yid);
            
            % Generate sample points within the population cell
            temp_lats = linspace(sample_lat - pop_lat_step / 2, sample_lat + pop_lat_step / 2, sample_step);
            temp_logs = linspace(sample_log - pop_log_step / 2, sample_log + pop_log_step / 2, sample_step);
            [mesh_lat, mesh_log] = meshgrid(temp_lats, temp_logs);
            
            % Check which of these sample points are inside the boundary
            % valid_points = inpolygon(mesh_log(:), mesh_lat(:), boundary.X, boundary.Y);
            
            % Loop through valid sample points and assign population to nearest nodes
            for idx = 1:length(mesh_lat(:))
                temp_lat = mesh_lat(idx);
                temp_log = mesh_log(idx);
                [node_id, ~] = identify_nearest_node_given_log_lat(node_data, temp_log, temp_lat);
                node_data(node_id, 4) = node_data(node_id, 4) + pop_Data(yid, xid) / (sample_step^2);
            end
        end
    end
end

% fprintf('\n');

% Ensure no negative population values in node_data
node_data(node_data(:, 4) < 0, 4) = 0;
node_data(:,1) = init_node_id;
% save(strcat(bd_path,system_type,'_node_data.mat'),'node_data');

% convhull(P(1).V)         
%     for r=1:length(raster_data(:,1))
%         rlog=raster_data(r,2);rlat=raster_data(r,3);
%         yid=floor((pop_R.LatitudeLimits(2)-rlat)/pop_lat_step)+1;xid=floor((rlog-pop_R.LongitudeLimits(1))/pop_log_step)+1;
%         
%         %generate a lot of population points around the raster area
%         sample_pop=zeros((2*ceil(lat_step/pop_lat_step)+1)*(2*ceil(log_step/pop_log_step)+1)*sample_step*sample_step,3);pk=0;
%         for klat=-ceil(lat_step/pop_lat_step):ceil(lat_step/pop_lat_step)
%             for klog=-ceil(log_step/pop_log_step):ceil(log_step/pop_log_step)
%                 lat_id=yid+klat;log_id=xid+klog;
%                 [sample_lat,sample_log]=intrinsicToGeographic(pop_R,log_id,lat_id);
%                 for i=1:sample_step
%                     temp_lat=sample_lat-pop_lat_step/2-pop_lat_step/(sample_step*2)+i*pop_lat_step/sample_step;
%                     for j=1:sample_step
%                         temp_log=sample_log-pop_log_step/2-pop_log_step/(sample_step*2)+j*pop_log_step/sample_step;
%                         pk=pk+1;sample_pop(pk,:)=[temp_lat temp_log pop_Data(lat_id,log_id)/(sample_step*sample_step)];
%                     end
%                 end
%             end
%         end
%         
%         
%         
%         
%         
%         raster_data(r,5)=sum(sample_pop(abs(sample_pop(:,1)-rlat)<=lat_step/2 & abs(sample_pop(:,2)-rlog)<=log_step/2,3));
%     end
% end
% 
% % figure;
% % for e=1:max(size(edge_str))
% %     ex=edge_str(e).X;ey=edge_str(e).Y;
% %     plot(ex,ey,'-','color',[0.5 0.5 0.5]);hold on;
% % end
end