%% After optimization
clear;
clc;
close all;

filename = input('Введите название файла без расширения: \n', 's');
basepath = 'C:\Users\User\Documents\MATLAB\AEP\Tasks and Scripts_v3\House_v2\Diplom\Workspaces\';
filepath = [basepath, filename, '.mat'];
load(filepath);
close all;
% Нужно два расклада по двум прогнозам:
%   Без исключения задач:
%     - Оптимистичный исход
%     - Пессимистичный исход
%   С исключением задач:
%     - Оптимистичный исход
%     - Пессимистичный исход


G_HPMa=digraph(PSM_hpma(:,1:N), Task_list_D);
G_HPMa_loops=digraph(PSM_hpma(:,1:N), Task_list_D, 'omitselfloop');

PSM_hpma_back=PSM_hpma(:,1:N);
D=diag(PSM_hpma_back);
% % D(:)=0;
% PSM_hpma_back=PSM_hpma_back+(diag(PSM_hpma_back)-D);

Start_points_HPMa = G_HPMa_loops.Nodes.Name(indegree(G_HPMa_loops,G_HPMa_loops.Nodes.Name)==0);
Finish_points_HPMa = G_HPMa_loops.Nodes.Name(outdegree(G_HPMa_loops,G_HPMa_loops.Nodes.Name)==0);

One_dots = string(intersect(Start_points_HPMa, Finish_points_HPMa));
Included_dots=Task_list_D(find(D));
Included_dots=Included_dots(ismember(Task_list_D(find(D)), One_dots));

% Start_points_HPMa = 
% Finish_points_HPMa = 

Start_tasks_op=string(G_HPMa_loops.Edges.EndNodes(:,1));
Next_tasks_op=string(G_HPMa_loops.Edges.EndNodes(:,2));

% Критические пути для обоих раскладов одинаковы отличаться будут лишь диаграммы ганта
[S_t_n_max, N_t_n_max, TrueTD_max, CritPath_op_max, CritTime_op_max] = True_Connects(Start_tasks_op, Next_tasks_op, Task_list_D, Tdmax);
[S_t_n_min, N_t_n_min, TrueTD_min, CritPath_op_min, CritTime_op_min] = True_Connects(Start_tasks_op, Next_tasks_op, Task_list_D, Tdmin);

G_HPMa_loops_nodots=G_HPMa_loops.rmnode(One_dots);

% Для логики добавим связь между первой точкой (началом всех работ) и критическим путём
if Start_tasks(1)~=CritPath_op_max(1) || Start_tasks(1)~=CritPath_op_min(1)

    % Для внутренних стартовых точек графа
    % Если есть "стартовые" задачи прилегающие к критическому пути -
    % создать для них связи от предыдущей задачи критического пути

    Start_tasks_op=cat(1,Start_points, Start_tasks_op);
    Next_tasks_op=cat(1,CritPath_op_max(1),Next_tasks_op);
        
    % Повторный поиск критического времени и пути на основе лишь одной стартовой задачи
    [S_t_n_max, N_t_n_max, TrueTD_max, CritPath_op_max, CritTime_op_max] = True_Connects(Start_tasks_op, Next_tasks_op, Task_list_D, Tdmax);
    [S_t_n_min, N_t_n_min, TrueTD_min, CritPath_op_min, CritTime_op_min] = True_Connects(Start_tasks_op, Next_tasks_op, Task_list_D, Tdmin);
    
    fprintf("Critical time max before adaptation: %3.2f\n", CritTime_op_max);
    fprintf("Critical time min before adaptation: %3.2f\n", CritTime_op_min);

    for i=2:length(CritPath_op_max)
        
        P=predecessors(G_HPMa_loops_nodots,string(CritPath_op_max(i)));
        P=P(ismember(P,Start_points_HPMa));
        if isempty(P) | ismember(P, CritPath_op_max)
            continue
        end
        Start_tasks_op=vertcat(Start_tasks_op, repmat(CritPath_op_max(i-1),length(P),1));
        Next_tasks_op=vertcat(Next_tasks_op, P);
        T=table(repmat(CritPath_op_max(i-1),length(P),1),P,ones(length(P),1));
        T.EndNodes=[T.Var1 T.P]; T.Weight=T.Var3;
        T(:, [1:3])=[];
        G_HPMa_loops_nodots=G_HPMa_loops_nodots.addedge(T);
        T=table([CritPath_op_max(1) CritPath_op_max(2)], 1,'VariableNames', {'EndNodes', 'Weight'});
        G_HPMa_loops_nodots=G_HPMa_loops_nodots.addedge(T);
       
    end
 
end

for i=1:length(Included_dots)
    % Происходит обход по критическому пути в обратном направлении
    count_down=length(CritPath_op_max)-1;

    % Пока длительность i-ой одиночной задачи не будет меньше
    % текущей задачи критического пути, то продолжаем обратный его обход
    while Tdmax(Task_list_D==Included_dots(i))>Tdmax(Task_list_D==string(CritPath_op_max(count_down)))
        count_down=count_down-1;
    end

    % И цепляем i-ую одиночную задачу к подходящей задаче на критическом пути
    % К векторам парных задач
    Start_tasks_op(end+1)=string(CritPath_op_max(count_down-1));
    Next_tasks_op(end+1)=Included_dots(i);

    % И к графу (переделать под любое количество!!)
    T.EndNodes(1)=cellstr(Start_tasks_op(end));
    T.EndNodes(2)=cellstr(Next_tasks_op(end));
    
    G_HPMa_loops_nodots=addedge(G_HPMa_loops_nodots,T);

end

% Повторный проход по критическому пути при одной стартовой и одной финишной задаче

Start_points_HPMa = G_HPMa_loops_nodots.Nodes.Name(indegree(G_HPMa_loops_nodots,G_HPMa_loops_nodots.Nodes.Name)==0);
Finish_points_HPMa = G_HPMa_loops_nodots.Nodes.Name(outdegree(G_HPMa_loops_nodots,G_HPMa_loops_nodots.Nodes.Name)==0);

T=[];
% Проход для "финишных" задач 
for i=1:length(Finish_points_HPMa)
    if ~isequal(Finish_points_HPMa(i),Finish_points)
        T=vertcat(T, table([Finish_points_HPMa(i) Finish_points], 1,'VariableNames', {'EndNodes', 'Weight'}));
        Start_tasks_op(end+1)=string(Finish_points_HPMa(i));
        Next_tasks_op(end+1)=string(Finish_points);
    end
end
% Проход для "стартовых" задач 
if numel(Start_points_HPMa)>1
    for i=1:length(Start_points_HPMa)
        if ~isequal(Start_points_HPMa(i),Start_points)
            T=vertcat(T, table([Start_points Start_points_HPMa(i)], 1,'VariableNames', {'EndNodes', 'Weight'}));
            Start_tasks_op(end+1)=string(Start_points);
            Next_tasks_op(end+1)=string(Start_points_HPMa(i));
        end
    end
end

T(1,:)=[];
G_HPMa_loops_nodots=addedge(G_HPMa_loops_nodots,T);

Start_points_HPMa = G_HPMa_loops_nodots.Nodes.Name(indegree(G_HPMa_loops_nodots,G_HPMa_loops_nodots.Nodes.Name)==0);
Finish_points_HPMa = G_HPMa_loops_nodots.Nodes.Name(outdegree(G_HPMa_loops_nodots,G_HPMa_loops_nodots.Nodes.Name)==0);

[S_t_n_max, N_t_n_max, TrueTD_max, CritPath_op_max, CritTime_op_max] = True_Connects(Start_tasks_op, Next_tasks_op, Task_list_D, Tdmax);
[S_t_n_min, N_t_n_min, TrueTD_min, CritPath_op_min, CritTime_op_min] = True_Connects(Start_tasks_op, Next_tasks_op, Task_list_D, Tdmin);

fprintf("Critical time max after adaptation: %3.2f\n", CritTime_op_max);
fprintf("Critical time min after adaptation: %3.2f\n", CritTime_op_min);

% Убираем из одиночных задач те, что включены в проект
% One_dots(One_dots==Included_dots)=[];
One_dots(ismember(One_dots,Included_dots))=[];

%% Dates

% Name=Task_list_D;
% GG=digraph(G_Vis.Edges, table(Name));

% Настройка списка задач для функции Timing (УМЕНЬШАЕТСЯ КОЛ-ВО ЗАДАЧ. ВЕРНУТЬ!)

Dels=ismember(Task_list_D, One_dots);

Task_list_D_op_nodots=Task_list_D;
Task_list_D_op_nodots(Dels)=[];
% return
% Тоже уменьшается!!
Tdmax_op=Tdmax;
Tdmin_op=Tdmin;
Tdmax_op(Dels)=[];
Tdmin_op(Dels)=[];

[S_t_n_max, N_t_n_max, TrueTD_max, CritPath_op_max, CritTime_op_max] = True_Connects(Start_tasks_op, Next_tasks_op, Task_list_D_op_nodots, Tdmax_op);
[S_t_n_min, N_t_n_min, TrueTD_min, CritPath_op_min, CritTime_op_min] = True_Connects(Start_tasks_op, Next_tasks_op, Task_list_D_op_nodots, Tdmin_op);

fprintf("Critical time max after adaptation: %3.2f\n", CritTime_op_max);
fprintf("Critical time min after adaptation: %3.2f\n", CritTime_op_min);
% return
%% Graphs

figure("Name",'Logic network of HPMa with excluded tasks','units','normalized','outerposition',[0 0 1 1])
ADD_Vis1=plot(G_HPMa, 'Layout','layered', 'AssignLayers','alap' , 'Direction', 'right');
title('Logic network of HPMa with excluded tasks');
ADD_Vis1.ArrowSize=6; % Increase arrows
ADD_Vis1.ArrowPosition=0.9;
%% ЗДЕСЬ ВЫЗВАТЬ ВИЗ
G_HPMa_loops=rmnode(G_HPMa_loops, One_dots);
% Name1='Project Graph of HPMa without excluded tasks (max TD)';
figure("Name",'Project Graph of HPMa without excluded tasks','units','normalized','outerposition',[0 0 1 1])
ADD_Vis2=plot(G_HPMa_loops, 'Layout','layered', 'AssignLayers','alap' , 'Direction', 'right');
title('Logic network of HPMa with excluded tasks');
ADD_Vis2.ArrowSize=6; % Increase arrows
ADD_Vis2.ArrowPosition=0.9;

% ADD_Vis(G_HPMa_loops, Start_tasks_op, Next_tasks_op, CritPath_op_min, Name1);
% ADD_Vis2_max=ADD_Vis(G_HPMa_loops, Start_tasks_op, Next_tasks_op, CritPath_op_max, Name1);

%% И ЗДЕСЬ

% figure("Name",'Adapted Project Graph of HPMa without excluded tasks (max TD)','units','normalized','outerposition',[0 0 1 1])
Name1='Adapted Project Graph of HPMa without excluded tasks (Td_{max})';
Name2='Adapted Project Graph of HPMa without excluded tasks (Td_{min})';

G_HPMa_loops_nodots=addedge(G_HPMa_loops_nodots, '2', '4', 1);
Finish_points_HPMa = {'24'};
ADD_Vis3_max=ADD_Vis(G_HPMa_loops_nodots, Start_points_HPMa, Finish_points_HPMa, CritPath_op_max, Name1, 'red');
ADD_Vis3_min=ADD_Vis(G_HPMa_loops_nodots, Start_points_HPMa, Finish_points_HPMa, CritPath_op_min, Name2, 'blue');
% Vis3=plot(G_HPMa_loops_nodots, 'Layout','layered', 'AssignLayers','alap' , 'Direction', 'right');
% title('Adapted Project Graph of HPMa without excluded tasks (max TD)');

%% Gantt chart (Charts)
Start_Date_Max_op = repmat(datetime(NaT), length(Task_list_D_op_nodots), 1);
Start_Date_Min_op = repmat(datetime(NaT), length(Task_list_D_op_nodots), 1);
Finish_Date_Max_op = repmat(datetime(NaT), length(Task_list_D_op_nodots), 1);
Finish_Date_Min_op = repmat(datetime(NaT), length(Task_list_D_op_nodots), 1);
Start_Date_Min_op(1)=datetime("01-Jan-2023");
Start_Date_Max_op(1)=Start_Date_Min_op(1);

[Start_Date_Max_op, Finish_Date_Max_op, Start_Date_Min_op, Finish_Date_Min_op] =...
    Timing(Start_Date_Max_op, Finish_Date_Max_op, Tdmax_op, Start_Date_Min_op, Finish_Date_Min_op, Tdmin_op, G_HPMa_loops_nodots, Task_list_D_op_nodots, CritPath_op_max, CritPath_op_min);

% Построение Диаграммы Ганта
figure("Name", 'Gantt Chart of HPMa (Td_{min})','units','normalized','outerposition',[0 0 1 1])
h = ganttChart(Task_list_D_op_nodots, Start_Date_Min_op, Finish_Date_Min_op);
h.TimeAxisLabel = 'Date';
h.TaskAxisLabel = 'Tasks';
h.Title = 'Gantt Chart of HPMa (Td_{min})';
h.FaceColor = "blue";
% h.Grid = 'on';

figure("Name",'Gantt Chart of HPMa (Td_{max})','units','normalized','outerposition',[0 0 1 1])
H = ganttChart(Task_list_D_op_nodots, Start_Date_Max_op, Finish_Date_Max_op);
H.TimeAxisLabel = 'Date';
H.TaskAxisLabel = 'Tasks';
H.Title = 'Gantt Chart of HPMa (Td_{max})';
H.FaceColor = "red";
% H.Grid = 'on';

lim_for_bars=ceil(CritTime/100)*100;
colors =['g', 'r'];

figure('Name','Durations','units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
B1=bar(categorical({'Optimistic','Pessimistic'}),[CritTime_f,CritTime]);
title('Duration before optimization')
ylim([0 lim_for_bars])

subplot(1,2,2)
B2=bar(categorical({'Optimistic','Pessimistic'}),[CritTime_op_min,CritTime_op_max]);
title('Duration after optimization')
ylim([0 lim_for_bars])

function Vis = ADD_Vis(G, Start_points, Finish_points, CritPath, Name, color)

figure("Name",Name,'units','normalized','outerposition',[0 0 1 1])
Vis=plot(G, 'Layout','layered', 'AssignLayers','alap' , 'Direction', 'right');
title(Name);
Vis.ArrowSize=6; % Increase arrows
highlight(Vis, G,'EdgeColor',color, 'NodeColor', 'red');
highlight(Vis,Start_points,"NodeColor","green",MarkerSize=4); % Highliting Start points
highlight(Vis,Finish_points,"NodeColor","cyan",MarkerSize=4); % Highliting Finish points
highlight(Vis, CritPath,"EdgeColor",'magenta', 'ArrowSize',10);
Vis.ArrowPosition=0.9;


end