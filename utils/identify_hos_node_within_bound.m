function hos_node = identify_hos_node_within_bound(hos_lon_lat,node_data,boundary)
    % IDENTIFY_HOS_NODE_WITHIN_BOUND identify the node_id of the nodes nearest to
    %    hosiptals within a given boundary. 
    % 
    % Example:
    % hos_node = identify_hos_node_within_bound(hos_lon_lat,node_data,boundary)
    %
    % Input: hos_lon_lat [lons, lats], the location of hospitals
    %        node_data [node_id, lon, lat], the nodes on the road network
    %        boundary boundary.X; boundary.Y; 
    % Output:
    %        hos_node [hospital_id, lon, lat, node_id] node_id: the id of 
    % road node corresponding to the hosiptal

    % load(strcat(hp_path,'china_hospital.mat'),'hospital');china_hos=hospital;
    % load(strcat(bd_path,num2str(city_id),city_name,'\sqr60_rs1_area_raster_info.mat'),'area_raster_info');
    
    
    % if strcmp(system_type,'admin')
    %     boundary=shaperead(strcat(osmnx_path,num2str(city_id),city_name,'\admin_boundary\boundary.shp'));
    % end
    % if strcmp(system_type,'sqr60')
    %     boundary.X=[clog-log_lat_step(1)*30.5 clog+log_lat_step(1)*30.5 clog+log_lat_step(1)*30.5 clog-log_lat_step(1)*30.5 NaN];
    %     boundary.Y=[clat-log_lat_step(2)*30.5 clat-log_lat_step(2)*30.5 clat+log_lat_step(2)*30.5 clat+log_lat_step(2)*30.5 NaN];
    % end

    if nargin < 3
        city_hos = hos_lon_lat;
    else
        city_hos = hos_lon_lat(inpolygon(hos_lon_lat(:,1),hos_lon_lat(:,2),boundary.X,boundary.Y)==1,:);
    end
    
    hos_node = zeros(length(city_hos(:,1)),1);
    for h=1:length(city_hos(:,1))
        [hos_node(h),~]=identify_nearest_node_given_log_lat(node_data,city_hos(h,1),city_hos(h,2));
    end
    
    hos_node = [[1:length(city_hos(:,1))]', city_hos, hos_node];
end