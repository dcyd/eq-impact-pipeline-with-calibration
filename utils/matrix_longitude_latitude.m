function distance=matrix_longitude_latitude(node_lon,node_lat,node_set_lon,node_set_lat)
% MATRIX_LONGITUDE_LATITUDE compute the distance from one node to a set of
% nodes
%
% distance=matrix_longitude_latitude(node_lon,node_lat,node_set_lon,node_set_lat)
%
% Input: [node_lon,node_at] [node_set_lon, node_set_lat]
% Output: distances from one node to node_set

node_lon=node_lon*3.14/180;
node_lat=node_lat*3.14/180;
node_set_lon=node_set_lon*3.14/180;
node_set_lat=node_set_lat*3.14/180;

temp=sqrt((sin((node_lat-node_set_lat)/2)).^2+cos(node_lat).*cos(node_set_lat).*(sin((node_lon-node_set_lon)/2)).^2);
distance=2*asin(temp)*6378.137;

