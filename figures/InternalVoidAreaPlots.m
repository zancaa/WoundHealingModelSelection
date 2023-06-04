%%
% This is a matlab script to import the void area and plot them.
% Note, must be in the InternalVoid/Circle results directory.

%%
% Import data
% 
% meshNoGhostInfRaw = importdata('Mesh/NoGhosts/InfiniteVT/voidArea.dat');
% meshNoGhostInfData = meshNoGhostInfRaw.data;

meshNoGhostFinRaw = importdata('Mesh/NoGhosts/FiniteVT/voidArea.dat');
meshNoGhostFinData = meshNoGhostFinRaw.data;

meshGhostRaw = importdata('Mesh/Ghosts/Void/voidArea.dat');
meshGhostData = meshGhostRaw.data;

%%%

VertexJaggedRaw = importdata('Vertex/Jagged/voidArea.dat');
VertexJaggedData = VertexJaggedRaw.data;

VertexSmoothRaw = importdata('Vertex/Smooth/voidArea.dat');
VertexSmoothData = VertexSmoothRaw.data;

VertexCurvedRaw = importdata('Vertex/Curved/voidArea.dat');
VertexCurvedData = VertexCurvedRaw.data;

%%
% Clean up and organise data

timeData = meshNoGhostFinData(:,1);

% meshNoGhostInfArea = meshNoGhostInfData(:,2);
meshNoGhostFinArea = meshNoGhostFinData(:,3);
meshGhostArea = meshGhostData(:,3);

NodeDefaultArea = NodeDefaultData(:,2);
NodeLargeArea = NodeLargeData(:,2);
NodeSmallArea = NodeSmallData(:,2);

VertexJaggedArea = VertexJaggedData(:,3);
VertexSmoothArea = VertexSmoothData(:,3);
VertexCurvedArea = VertexCurvedData(:,3);

%%
% Plot the void area

mLineWidth = 2.5;

figure
hold on;
plot(VertexSmoothData(:,1),VertexSmoothData(:,3),'-','linewidth',2,'Color',(1/255)*[23 63 95])
plot(VertexJaggedData(1:10:end,1),VertexJaggedData(1:10:end,3),'--','linewidth',2,'Color',(1/255)*[23 63 95])
plot(VertexCurvedData(:,1),VertexCurvedData(:,3),':','linewidth',2,'Color',(1/255)*[23 63 95])
legend('Smooth','Default','Curved','fontsize',12,'interpreter','Latex',...
    'EdgeColor','None','Color','None')
xlabel('Time (hours)','fontsize',14,'interpreter','Latex')
ylabel('Wound Area (CD$^2$)','fontsize',14,'interpreter','Latex')
title('Vertex model','fontsize',16,'interpreter','Latex')

figure
hold on
plot(meshNoGhostFinData(:,1),zeros(length(meshNoGhostFinData(:,1)),1),'-','LineWidth',2,'Color',(1/255)*[237 85 59])
plot(meshNoGhostFinData(:,1),meshNoGhostFinData(:,3),'--','LineWidth',2,'Color',(1/255)*[237 85 59])
plot(meshGhostData(:,1),meshGhostData(:,3),':','LineWidth',2,'Color',(1/255)*[237 85 59])
% plot(timeData,meshNoGhostInfArea,'--','linewidth',1.5)
legend('Unbounded','Bounded','Ghost','fontsize',12,'interpreter','Latex',...
    'EdgeColor','None','Color','None')

xlabel('Time (hours)','fontsize',14,'interpreter','Latex')
ylabel('Wound Area (CD$^2$)','fontsize',14,'interpreter','Latex')
title('Voronoi tessellation','fontsize',16,'interpreter','Latex')
% title('Void Area with Domain scaling of 0.8','fontsize',16,'interpreter','Latex')

%set(gca, 'xscale','log')

axis([0.49 timeData(end) -0.01 25])

%% Calculate wound area for node simulation

% Repulsion
%%%
NodeSmallRaw = importdata('Node/SmallCutoff/Post-Void/voidArea.dat');
NodeSmallData = NodeSmallRaw.data;
NodeSmallData = [NodeSmallData zeros(length(NodeSmallData(:,1)),1)];
for i = 1:length(NodeSmallData(:,1))
    filename1 = 'Node/SmallCutoff/Masks/Repulsion_mask_t_';
    if i < 11
        filename = strcat(filename1,'00',num2str(i-1),'.png');
    elseif i < 101
        filename = strcat(filename1,'0',num2str(i-1),'.png');
    else
        filename = strcat(filename1,num2str(i-1),'.png');
    end
    BW = imread(filename);
    I = rgb2gray(BW);
    area_tmp = bwarea(I); %pixel area
    area_perc = area_tmp/(536*644); %percentage area 
    NodeSmallData(i,4) = area_perc*16*16*sqrt(3)*0.5; %multiply percentage by domain size
end


% Default
NodeDefaultRaw = importdata('Node/DefaultCutoff/Post-Void/voidArea.dat');
NodeDefaultData = NodeDefaultRaw.data;
NodeDefaultData = [NodeDefaultData zeros(length(NodeDefaultData(:,1)),1)];
for i = 1:length(NodeDefaultData(:,1))
    filename1 = 'Node/DefaultCutOff/Masks/Short_range_mask_t_';
    if i < 11
        filename = strcat(filename1,'00',num2str(i-1),'.png');
    elseif i < 101
        filename = strcat(filename1,'0',num2str(i-1),'.png');
    else
        filename = strcat(filename1,num2str(i-1),'.png');
    end
    BW = imread(filename);
    I = rgb2gray(BW);
    area_tmp = bwarea(I); %pixel area
    area_perc = area_tmp/(536*644); %percentage area 
    NodeDefaultData(i,4) = area_perc*16*16*sqrt(3)*0.5; %multiply percentage by domain size
end

% Long-range
NodeLargeRaw = importdata('Node/LargeCutoff/Post-Void/voidArea.dat');
NodeLargeData = NodeLargeRaw.data;
NodeLargeData = [NodeLargeData zeros(length(NodeLargeData(:,1)),1)];
for i = 1:length(NodeLargeData(:,1))
    filename1 = 'Node/LargeCutoff/Masks/Long_range_mask_t_';
    if i < 11
        filename = strcat(filename1,'00',num2str(i-1),'.png');
    elseif i < 101
        filename = strcat(filename1,'0',num2str(i-1),'.png');
    else
        filename = strcat(filename1,num2str(i-1),'.png');
    end
    BW = imread(filename);
    I = rgb2gray(BW);
    area_tmp = bwarea(I); %pixel area
    area_perc = area_tmp/(536*644); %percentage area 
    NodeLargeData(i,4) = area_perc*16*16*sqrt(3)*0.5; %multiply percentage by domain size
end

figure
hold on
plot(NodeSmallData(:,1),NodeSmallData(:,4),'-','linewidth',2,'Color',(1/255)*[60 174 163])
plot(NodeDefaultData(:,1),NodeDefaultData(:,4),'--','linewidth',2,'Color',(1/255)*[60 174 163])
plot(NodeLargeData(:,1),NodeLargeData(:,4),':','linewidth',2,'Color',(1/255)*[60 174 163])
legend('Repulsion','Short-range','Long-range','fontsize',12,'interpreter','Latex',...
    'EdgeColor','None','Color','None')
xlabel('Time (hours)','fontsize',14,'interpreter','Latex')
ylabel('Wound Area (CD$^2$)','fontsize',14,'interpreter','Latex')
title('Overlapping spheres','fontsize',16,'interpreter','Latex')