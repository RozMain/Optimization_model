function [] = Visualising(G_Vis, g_Vis,Start_points, Finish_points, TaskType, GName, GName2, CritTime_max, CritTime_min, CritPath, CritPath_f)
%% Отрисовка графа с учётом гибкости
figure('Name', 'Graph with Td_{max}','units','normalized','outerposition',[0 0 1 1])
Visual=plot(G_Vis,'Layout','layered', 'AssignLayers','alap' , 'Direction', 'right');
title(GName);
Visual.ArrowSize=4; % Increase arrows
highlight(Visual, G_Vis,'EdgeColor','red', 'NodeColor', 'red');
highlight(Visual,Start_points,"NodeColor","green",MarkerSize=4); % Highliting Start points
highlight(Visual,Finish_points,"NodeColor","cyan",MarkerSize=4); % Highliting Finish points
highlight(Visual, CritPath,"EdgeColor",'magenta', 'ArrowSize',10);

%% Отрисовка графа обязательных задач и зависимостей
figure ('Name', 'Graph with Td_{min}','units','normalized','outerposition',[0 0 1 1])
Visual2=plot(g_Vis,'Layout','layered', 'AssignLayers','alap', 'Direction', 'right');
Visual2.ArrowSize=4;
title(GName2);
highlight(Visual2,Start_points,"NodeColor","green",MarkerSize=4); % Highliting Start points
highlight(Visual2,Finish_points,"NodeColor","cyan",MarkerSize=4); % Highliting Finish points
highlight(Visual2, CritPath_f,"EdgeColor",'magenta', 'ArrowSize',10);

figure
bar(categorical({'TPTmin', 'TPTmax'}),[CritTime_min, CritTime_max])
title('Duration')

%% Отрисовка графа с классификацией задач

figure ('Name', 'Classification','units','normalized','outerposition',[0 0 1 1])
hold on
scatter(0,0.1,"black", "filled");
scatter(0,0.1,"red", "filled");
scatter(0,0.1,"blue", "filled");


Visual=plot(g_Vis,'Layout','layered', 'AssignLayers','alap' , 'Direction', 'right');
title('House (Classification of tasks)');
Visual.ArrowSize=4; % Increase arrows
highlight(Visual, g_Vis,'EdgeColor','red', 'NodeColor', 'red');
% highlight(Visual, CritPath,"EdgeColor",'magenta', 'ArrowSize',10);
hold off

B=[];
R=[];
Bl=[];
Y=[];

for i=1:length(TaskType)
    switch TaskType(i)
        case 11
            B = vertcat(B, i);
        case 22
            R = vertcat(R, i);
        case 33
            Bl = vertcat(Bl, i);
        case 44
            Y = vertcat(Y, i);
    end
end
legend({'Строительсто и монтажные работы', 'Работы обеспечивающие коммуникации', 'Работа с документацией', 'Работы связанные с настройкой, диагностикой и пуско-наладкой', 'Граф'}, 'Location','southwest');
highlight(Visual,B,"NodeColor","black");
% legend('Строительсто и монтажные работы')
highlight(Visual,R,"NodeColor","red");
% legend('Работы обеспечивающие коммуникации')
highlight(Visual,Bl,"NodeColor","blue");
% legend('Работа с документацией');
highlight(Visual,Y,"NodeColor","yellow");
% legend('Работы связанные с настройкой, диагностикой и пуско-наладкой');
clear GName
clear GName2
end
