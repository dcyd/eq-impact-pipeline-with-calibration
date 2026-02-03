function [node_id,edist]=identify_nearest_node_given_log_lat(node_data,given_log,given_lat)
%IDENTIFY_NEAREST_GIVEN_LOG_LAT identify the nearest node for a given point
%among the node data 
% EXAMPLE:
%   [node_id,edist]=identify_nearest_node_given_log_lat(node_data,given_log,given_lat)
% INPUT
%     node_data
%     given_log
%     given_lat
% OUTPUT
%   node_id
%   edist (km)

define_constants;
kcx = given_log;
kcy = given_lat;
mdist = min(abs(node_data(:,NX)-kcx)+abs(node_data(:,NY)-kcy));
kk_node = node_data(abs(node_data(:,NX)-kcx)<=mdist*1.5 & abs(node_data(:,NY)-kcy)<=mdist*1.5,[NI NX NY NI]);
% for j=1:length(kk_node(:,1))
%     kk_node(j,4)=longitude_latitude(kcx,kcy,kk_node(j,NX),kk_node(j,NY));
% end
kk_node(:,4) = matrix_longitude_latitude(kcx,kcy,kk_node(:,NX),kk_node(:,NY));
kk_node = sortrows(kk_node,4);
node_id = kk_node(1,1);
edist = kk_node(1,4);
        
    