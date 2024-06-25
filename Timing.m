function [Start_Date_Max, Finish_Date_Max, Start_Date_Min, Finish_Date_Min] =...
    Timing(Start_Date_Max, Finish_Date_Max, Tdmax, Start_Date_Min, Finish_Date_Min, Tdmin, GG, Task_list_D, CritPath, CritPath_f)

for i=1:length(CritPath)

    Finish_Date_Max(find(CritPath{i}==Task_list_D))=Start_Date_Max(find(CritPath{i}==Task_list_D))+days(Tdmax(find(CritPath{i}==Task_list_D)));
    if i<length(CritPath)

        Start_Date_Max(find(CritPath{i+1}==Task_list_D))=Finish_Date_Max(find(CritPath{i}==Task_list_D));

    end
end

[Start_Date_Max, Finish_Date_Max] = FillDates_NearCritPath(Start_Date_Max, Finish_Date_Max, CritPath, Task_list_D, Tdmax, GG);

for i=1:length(CritPath_f)
    Finish_Date_Min(find(CritPath_f{i}==Task_list_D))=Start_Date_Min(find(CritPath_f{i}==Task_list_D))+days(Tdmin(find(CritPath_f{i}==Task_list_D)));
    if i<length(CritPath_f)
      
        Start_Date_Min(find(CritPath_f{i+1}==Task_list_D))=Finish_Date_Min(find(CritPath_f{i}==Task_list_D));

    end
end

[Start_Date_Min, Finish_Date_Min] = FillDates_NearCritPath(Start_Date_Min, Finish_Date_Min, CritPath_f, Task_list_D, Tdmin, GG);

    function [Start_Date, Finish_Date] = FillDates_NearCritPath(Start_Date, Finish_Date, Crit_Path, Task_list_D, Td, G)
    %% 2й проход для заполнения дат у задач прилегающих к критическому пути задачам
    for Crit_Task = Crit_Path
        % Найти последователей задачи
        NTasks = GetSuccessors(Task_list_D, Crit_Task, G);

        % Заполнить даты для каждого последователя
        for j = 1:length(NTasks)
            if any(Task_list_D(NTasks(j)) == Crit_Path)
                continue
            end

            Start_Date(NTasks(j)) = Finish_Date(Task_list_D == Crit_Task);
            Finish_Date(NTasks(j)) = Start_Date(NTasks(j)) + Td(NTasks(j));
        end
    end

    %% 3й проход для заполнения дат у оставшихся веток (возможно возникновение рекурсии!!!)
%     NTasks=[];

num=1;
    while any(isnat(Start_Date))
        n = find(isnat(Start_Date),num);
        if ~ismember(Task_list_D(n), string(table2array(G.Nodes)))
            continue
        end
        
        while isnat(Finish_Date(n))
            n=GetPredessors(Task_list_D, Task_list_D(n), G);
            NTasks = GetSuccessors(Task_list_D, Task_list_D(n), G);
            num=num+1;
        end
        num=1;

        for j = 1:length(NTasks)
            if ~isnat(Start_Date(NTasks(j)))
                continue
            end

            Start_Date(NTasks(j)) = Finish_Date(n);
            Finish_Date(NTasks(j)) = Start_Date(NTasks(j)) + Td(NTasks(j));
        end
    end
    end

    function NTasks = GetSuccessors(Task_list_D, Task, G)
    % Найти последователей задачи
    NTasks = [];
    successors = G.successors(Task);
    for j = 1:length(successors)
        NTasks = [NTasks; find(Task_list_D == successors(j))];
    end
    end

    function NTasks = GetPredessors(Task_list_D, Task, G)
    % Найти последователей задачи
    NTasks = [];
    successors = G.predecessors(Task);
    for j = 1:length(successors)
        NTasks = [NTasks; find(Task_list_D == successors(j))];
    end
    end

%     function [Start_Date, Finish_Date] = FillDates_NearCritPath(Start_Date, Finish_Date, Crit_Path, Task_list_D, Td, G)
%         %% 2й проход для заполнения дат у задач прилегающих к критическому пути задачам
%         for x=1:length(Crit_Path)
%             % Сколько последователей у задачи?
%             Tasks=G.successors(Task_list_D(Task_list_D==Crit_Path(x)));
%             NTasks=[];
%             for j=1:length(Tasks)
%                 NTasks=[NTasks; find(Task_list_D==Tasks(j))];
%             end
% 
%             % Универсальное решение для люого кол-ва последователей
%             if length(NTasks)>1
%                 for j=1:length(NTasks)
%                     if sum(Task_list_D(NTasks(j))==Crit_Path)
%                         continue
%                     end
% %                     Task_suc=Task_list_D(Task_list_D==G.successors(Task_list_D(NTasks(j))));
%                     Start_Date(NTasks(j))=Finish_Date(Task_list_D==Crit_Path(x));
%                     Finish_Date(NTasks(j))=Start_Date(NTasks(j))+Td(NTasks(j));
%                 end
%             end
%         end
%         
%         %% 3й проход для заполнения дат у оставшихся веток (возможно возникновение рекурсии!!!) 
%         % НАЙТИ ТОЧКУ ВОЗНИКНОВЕНИЯ РЕКУРСИИ И ИСПОЛЬЗОВАТЬ WHILE
% 
%         n=0;
%         while sum(isnat(Start_Date))~=0
%             n=n+1;
%             if ~ismember(Task_list_D(n), string(table2array(G.Nodes)))
%                 continue
%             end
%             Tasks=G.successors(Task_list_D(n));
%             NTasks=[];
% 
%             for j=1:length(Tasks)
%                 NTasks=[NTasks; find(Task_list_D==Tasks(j))];
%             end
%             
%             for j=1:length(NTasks)
%                 if ~isnat(Start_Date(NTasks(j)))
%                     continue
%                 end
%                 Start_Date(NTasks(j))=Finish_Date(n);
%                 Finish_Date(NTasks(j))=Start_Date(NTasks(j))+Td(NTasks(j));
%             end
% 
%             if n==length(Start_Date)
%                 n=0;
%             end
%         end
%     end

end
% end