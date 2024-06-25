function [Start_tasks, Next_tasks, ProbabilityT, ProbabilityC, Task_list_D, Tdmin, Tdmax, Cdmin, Cdmax, Rdmin, Rdmax, TaskType] =...
    Reconfig(Start_tasks, Next_tasks, GG, ProbabilityT, ProbabilityC, Task_list_D, Tdmin, Tdmax, Cdmin, Cdmax, Rdmin, Rdmax, TaskType)

codes = {'a' 'b' 'c'};


% Q0="Сколько у вас задач, имеющих альтернативы (помеченные буквами a,b и c)?\n";
% A0=input(Q0);

NtasksA=find(contains(Task_list_D,'a'));
tasksA=Task_list_D(NtasksA);
NtasksB=find(contains(Task_list_D,'b'));
tasksB=Task_list_D(NtasksB);
NtasksC=find(contains(Task_list_D,'c'));
tasksC=Task_list_D(NtasksC);


Rele=[tasksA'; tasksB'; tasksC'];
Rele=allcombinations(Rele);
% r - кол-во комбинаций r==(length(codes)^length(tasksA));
% с - кол-во задач с альтернативами
[r c]=size(Rele);

% n_chk_s=strip(tasksA,"right","a");


%Дореализовать случаи если не будет одной из букв!!
% if length(tasksA)==length(tasksB) && length(tasksA)==length(tasksC)
%     counters(1:length(tasksA),1)='a';
% else
%     [~,N]=max([length(tasksA), length(tasksB), length(tasksC)]);
%     
%     switch N
%         case 1
%             
%         case 2
% 
%         case 3
% 
%     end
%     
% end

n=0;
% СТРЕМИМСЯ К WHILE ДЛЯ РАСШИРЯЕМОСТИ!!!
while ~(n==r)
    n=n+1;
    
    task_for_del=[];
    N_del=[];
    name=[];

    name=Rele(n,:);

    %     name=[name Rele{n} Rele{r+n}];
    %     name=name';
    %     chk_symb=erase(name,n_chk_s);

    task_for_del=setdiff([tasksA tasksB tasksC], string(Rele(n,:)));
    N_del=find(contains(Task_list_D, task_for_del));

    %     for i=1:c
%
%                switch chk_symb{i}
%             case 'a'
%                 task_for_del=[task_for_del; tasksB(i)];
%                 task_for_del=[task_for_del; tasksC(i)];
%
%                 N_del=[N_del; NtasksB(i)];
%                 N_del=[N_del; NtasksC(i)];
%
%             case 'b'
%                 task_for_del=[task_for_del; tasksA(i)];
%                 task_for_del=[task_for_del; tasksC(i)];
%
%                 N_del=[N_del; NtasksA(i)];
%                 N_del=[N_del; NtasksC(i)];
%
%             case 'c'
%                 task_for_del=[task_for_del; tasksA(i)];
%                 task_for_del=[task_for_del; tasksB(i)];
%
%                 N_del=[N_del; NtasksA(i)];
%                 N_del=[N_del; NtasksB(i)];
%
%         end
%     end


    [GG_mod,Tdmax_mod, Tdmin_mod, Cdmin_mod, Cdmax_mod, Rdmin_mod, Rdmax_mod, TaskType_mod,  ProbabilityT_mod, ProbabilityC_mod, S_t_mod, N_t_mod]=...
        DelNode(task_for_del, N_del, GG, Tdmax, Tdmin, Task_list_D, Cdmin, Cdmax, Rdmin, Rdmax, TaskType,  ProbabilityT, ProbabilityC, Start_tasks, Next_tasks);
    Tasks=string((GG_mod.Edges{:,:}));
    TlD_mod=string(table2array(GG_mod.Nodes));

    [~, ~, ~, CritPath_max, CritTime_max] = True_Connects(Tasks(1:end,1) , Tasks(1:end,2), TlD_mod, Tdmax_mod);
    [~, ~, ~, CritPath_min, CritTime_min] = True_Connects(Tasks(1:end,1) , Tasks(1:end,2), TlD_mod, Tdmin_mod);
    aaa=1;
    Struct_G(n,:)={string([name{1,:}]),GG_mod, CritTime_min, CritTime_max, CritPath_min, CritPath_max,...
        Tdmin_mod, Tdmax_mod, Cdmin_mod, Cdmax_mod, Rdmin_mod, Rdmax_mod, TaskType_mod,  ProbabilityT_mod, ProbabilityC_mod, TlD_mod, S_t_mod, N_t_mod,...
        sum(Cdmin_mod), sum(Cdmax_mod), sum(Rdmin_mod), sum(Rdmax_mod)};

end
A_V=2;
if c>2
    Q_V="\nВнимание! В вашем проекте много альтернативных вариантов.\n" + ...
        "Произвести визуализацию всех сценариев? (Может потребоваться какое то время и графы будут практически не читабельными!):\n" + ...
        "1 - Нет\n" + ...
        "2 - Всё равно произвести визуализацию\n";
    A_V=input(Q_V);
end

if A_V==2
    figure('Name','All project scenarios','units','normalized','outerposition',[0 0 1 1]);
    for i=1:length(codes)^c
        subplot((length(codes)^c)/length(codes),length(codes),i)
        Visual=plot(Struct_G{i,2}, 'Layout','layered', 'AssignLayers','alap' , 'Direction', 'right');
        title(Struct_G{i,1});
        txtfs=['Critical time min = {\color{green}' num2str(Struct_G{i,3})  '} days. Critical time max = {\color{red}' num2str(Struct_G{i,4}) '} days'];
        subtitle(txtfs, "BackgroundColor","white", "FontWeight","bold")

        highlight(Visual, Struct_G{i,5},"EdgeColor",'green', 'ArrowSize',10);
        highlight(Visual, Struct_G{i,6},"EdgeColor",'red', 'ArrowSize',10);
    end
end

Q1_1="\nВпишите ограничение по бюджету (в тыс руб):\n";
CostCup=input(Q1_1);

% Sum_Cost=sum(Struct_G{:,10});
Q1="\nПровести предварительный поиск оптимального решения?\n" + ...
    "1- Да\n" + ...
    "2- Нет\n";
A1=input(Q1);
while isstring(A1) || ischar(A1) || A1 > 2 || A1 < 1
    if isstring(A1) || ischar(A1)
        fprintf("\nВведён недопустимый знак, повторите попытку!\n")
        A1=input(Q1);
    elseif A1 > 2 || A1 < 1
        fprintf("\nВведено некорректное значение, повторите попытку!\n")
        A1=input(Q1);
    end
end


while isstring(CostCup) || ischar(CostCup) || CostCup<0
    if isstring(CostCup) || ischar(CostCup)
        fprintf("\nВведён недопустимый знак! Введите число!\n")
        CostCup=input(Q1_1);
    elseif CostCup<0
        fprintf("\nВведено отрицательное значение! Повторите попытку!\n")
        CostCup=input(Q1_1);
    end
end

if A1==2
    Q2="\nВыберите стратегию планирования:\n" + ...
        "1 - Уменьшение длительности проекта (может привести к увеличению стоимости)\n" + ...
        "2 - Уменьшение стоимости проекта (может привести к увеличению длительности)\n" + ...
        "3 - Уменьшение ресурсных затрат (может привести к увеличению стоимости)\n" + ...
        "4 - Оптимальное решение\n";
    A2 = input(Q2);

    while isstring(A2) || ischar(A2) || A2 > 4 || A2 < 1

        if isstring(A2) || ischar(A2)
            fprintf("\nВведён недопустимый знак, повторите попытку!\n")
            A2 = input(Q2);
        elseif A2 > 4 || A2 < 1
            fprintf("\nВведено некорректное значение, повторите попытку!\n")
            A2 = input(Q2);
        end
    end

    New_Struct_G = Cases(Struct_G,A2, CostCup);
    
    
    if isempty(New_Struct_G) && A2==4
        fprintf("\nОптимального автономного решения не выявлено! " + ...
            "Воспользуйтесь одним из частных решений.\n")

        Q2="\nВыберите стратегию планирования:\n" + ...
            "1 - Уменьшение длительности проекта (может привести к увеличению стоимости)\n" + ...
            "2 - Уменьшение стоимости проекта (может привести к увеличению длительности)\n" + ...
            "3 - Уменьшение ресурсных затрат (может привести к увеличению стоимости)\n";
        A2 = input(Q2);

        while isstring(A2) || ischar(A2) || A2 > 3 || A2 < 1

            if isstring(A2) || ischar(A2)
                fprintf("\nВведён недопустимый знак, повторите попытку!\n")
                A2 = input(Q2);
            elseif A2 > 3 || A2 < 1
                fprintf("\nВведено некорректное значение, повторите попытку!\n")
                A2 = input(Q2);
            end
        end
        New_Struct_G = Cases(Struct_G,A2, CostCup);
     end
else
    
    New_Struct_G = Cases(Struct_G,4, CostCup);
    if isempty(New_Struct_G)
        fprintf("\nОптимального автономного решения не выявлено! " + ...
            "Воспользуйтесь одним из частных решений.\n")
    end

    Q2="\nВыберите стратегию планирования:\n" + ...
        "1 - Уменьшение длительности проекта (может привести к увеличению стоимости)\n" + ...
        "2 - Уменьшение стоимости проекта (может привести к увеличению длительности)\n" + ...
        "3 - Уменьшение ресурсных затрат (может привести к увеличению стоимости)\n";
    A2 = input(Q2);

    while isstring(A2) || ischar(A2) || A2 > 3 || A2 < 1

        if isstring(A2) || ischar(A2)
            fprintf("\nВведён недопустимый знак, повторите попытку!\n")
            A2 = input(Q2);
        elseif A2 > 3 || A2 < 1
            fprintf("\nВведено некорректное значение, повторите попытку!\n")
            A2 = input(Q2);
        end
    end
    New_Struct_G = Cases(Struct_G,A2, CostCup);
end



% New_Struct_G = Cases(Struct_G,A2, CostCup);

[Start_tasks, Next_tasks, ProbabilityT, ProbabilityC, Task_list_D, Tdmin, Tdmax, Cdmin, Cdmax, Rdmin, Rdmax, TaskType] =...
                    Filling_exit(New_Struct_G);


%Доделать!!!
   
    function New_Struct_G = Cases(Struct_G,check, CostCup)

        switch check
            % Довести до ума алгоритм выбора!!!
            case 1

                Number_of_Strategy=[Struct_G{:,20}]<=CostCup; % Исключение вариантов превышения бюджета
                New_Struct_G=Struct_G(Number_of_Strategy,:);

                Number_of_Strategy=[New_Struct_G{:,4}]==min([New_Struct_G{:,4}]); % Фильтр по CritTime_max
                New_Struct_G=New_Struct_G(Number_of_Strategy,:);

                if sum(Number_of_Strategy)>1
                    Number_of_Strategy=[New_Struct_G{:,3}]==min([New_Struct_G{:,3}]); % Фильтр по CritTime_min
                    New_Struct_G=New_Struct_G(Number_of_Strategy,:);
                    if sum(Number_of_Strategy)>1
                        Number_of_Strategy=[New_Struct_G{:,19}]==max([New_Struct_G{:,19}]); % Фильтр по Cdmin
                        New_Struct_G=New_Struct_G(Number_of_Strategy,:);
                        if sum(Number_of_Strategy)>1
                            Number_of_Strategy=[New_Struct_G{:,20}]==max([New_Struct_G{:,20}]); % Фильтр по Cdmax
                            New_Struct_G=New_Struct_G(Number_of_Strategy,:);
                        end
                    end
                end

            case 2

                Number_of_Strategy=[Struct_G{:,20}]<=CostCup; % Исключение вариантов превышения бюджета
                New_Struct_G=Struct_G(Number_of_Strategy,:);

                Number_of_Strategy=[New_Struct_G{:,20}]==min([New_Struct_G{:,20}]); % Фильтр по Cdmax
                New_Struct_G=New_Struct_G(Number_of_Strategy,:);

                if sum(Number_of_Strategy)>1
                    Number_of_Strategy=[New_Struct_G{:,19}]==min([New_Struct_G{:,19}]); % Фильтр по Cdmin
                    New_Struct_G=New_Struct_G(Number_of_Strategy,:);
                    if sum(Number_of_Strategy)>1
                        Number_of_Strategy=[New_Struct_G{:,3}]==max([New_Struct_G{:,3}]); % Фильтр по CritTime_min
                        New_Struct_G=New_Struct_G(Number_of_Strategy,:);
                        if sum(Number_of_Strategy)>1
                            Number_of_Strategy=[New_Struct_G{:,4}]==max([New_Struct_G{:,4}]); % Фильтр по CritTime_max
                            New_Struct_G=New_Struct_G(Number_of_Strategy,:);
                        end
                    end
                end

            case 3

                Number_of_Strategy=[Struct_G{:,20}]<=CostCup; % Исключение вариантов превышения бюджета
                New_Struct_G=Struct_G(Number_of_Strategy,:);

                Number_of_Strategy=[New_Struct_G{:,22}]==min([New_Struct_G{:,22}]); % Фильтр по Rdmax
                New_Struct_G=New_Struct_G(Number_of_Strategy,:);

                if sum(Number_of_Strategy)>1
                    Number_of_Strategy=[New_Struct_G{:,21}]==min([New_Struct_G{:,21}]); % Фильтр по Rdmin
                    New_Struct_G=New_Struct_G(Number_of_Strategy,:);
                    if sum(Number_of_Strategy)>1
                        Number_of_Strategy=[New_Struct_G{:,19}]==max([New_Struct_G{:,19}]); % Фильтр по Cdmin
                        New_Struct_G=New_Struct_G(Number_of_Strategy,:);
                        if sum(Number_of_Strategy)>1
                            Number_of_Strategy=[New_Struct_G{:,20}]==max([New_Struct_G{:,20}]); % Фильтр по Cdmax
                            New_Struct_G=New_Struct_G(Number_of_Strategy,:);
                        end
                    end
                end

            case 4

                Number_of_Strategy=[Struct_G{:,20}]<=CostCup; % Исключение вариантов превышения бюджета
                New_Struct_G=Struct_G(Number_of_Strategy,:);

                Number_of_Strategy_1=[New_Struct_G{:,4}]==min([New_Struct_G{:,4}]); % Фильтр по CritTime_max
                Number_of_Strategy_2=[New_Struct_G{:,20}]==min([New_Struct_G{:,20}]); % Фильтр Cdmax
                Number_of_Strategy_3=[New_Struct_G{:,22}]==min([New_Struct_G{:,22}]); % Фильтр Rdmax
                
                
                if sum(Number_of_Strategy_1.*Number_of_Strategy_2.*Number_of_Strategy_3)==0 &&...
                        sum(Number_of_Strategy_1.*Number_of_Strategy_2)==0 &&...
                        sum(Number_of_Strategy_1.*Number_of_Strategy_3)==0 &&...
                        sum(Number_of_Strategy_2.*Number_of_Strategy_3)==0
                    New_Struct_G=[];
                    return
                end
                if sum(Number_of_Strategy_1.*Number_of_Strategy_2.*Number_of_Strategy_3)==0
                    if sum(Number_of_Strategy_1.*Number_of_Strategy_2)==0
                        if sum(Number_of_Strategy_1.*Number_of_Strategy_3)==0
                            % Варианта: if
                            % sum(Number_of_Strategy_2.*Number_of_Strategy_3)==0 не
                            % нужно, так как это последнее вложение фильтрации
                            % и оно представлено ниже

                            Number_of_Strategy=logical(Number_of_Strategy_2.*Number_of_Strategy_3); % вот оно
                            New_Struct_G=Struct_G(Number_of_Strategy,:);
                            if sum(Number_of_Strategy)>1
                                Number_of_Strategy=[New_Struct_G{:,4}]==min([New_Struct_G{:,4}]);
                                New_Struct_G=New_Struct_G(Number_of_Strategy,:);
                            end


                        else
                            Number_of_Strategy=logical(Number_of_Strategy_1.*Number_of_Strategy_3);
                            New_Struct_G=Struct_G(Number_of_Strategy,:);
                            if sum(Number_of_Strategy)>1
                                Number_of_Strategy=[New_Struct_G{:,20}]==min([New_Struct_G{:,20}]);
                                New_Struct_G=New_Struct_G(Number_of_Strategy,:);
                            end

                        end
                    else
                        Number_of_Strategy=logical(Number_of_Strategy_1.*Number_of_Strategy_2);
                        New_Struct_G=Struct_G(Number_of_Strategy,:);
                        if sum(Number_of_Strategy)>1
                            Number_of_Strategy=[Struct_G{:,22}]==min([Struct_G{:,22}]);
                            New_Struct_G=New_Struct_G(Number_of_Strategy,:);
                        end

                    end
                else
                    Number_of_Strategy=logical(Number_of_Strategy_1.*Number_of_Strategy_2.*Number_of_Strategy_3);
                    New_Struct_G=Struct_G(Number_of_Strategy,:);
                    % Маловероятно, что одновременно по всем трём фильтрам
                    % получится больше одного решения
                end



        end

    end

    function [Start_tasks, Next_tasks, ProbabilityT, ProbabilityC, Task_list_D, Tdmin, Tdmax, Cdmin, Cdmax, Rdmin, Rdmax, TaskType]=...
            Filling_exit(Struct_G)

        Start_tasks=Struct_G{17};
        Next_tasks=Struct_G{18};
        ProbabilityT=Struct_G{14};
        ProbabilityC=Struct_G{15};
        Task_list_D=Struct_G{16};
        Tdmin=Struct_G{7};
        Tdmax=Struct_G{8};
        Cdmin=Struct_G{9};
        Cdmax=Struct_G{10};
        Rdmin=Struct_G{11};
        Rdmax=Struct_G{12};
        TaskType=Struct_G{13};

    end

    function [GG_mod,Tdmax_mod, Tdmin_mod, Cdmin_mod, Cdmax_mod, Rdmin_mod, Rdmax_mod, TaskType_mod,  ProbabilityT_mod, ProbabilityC_mod, S_t_mod, N_t_mod] =...
            DelNode(task_for_del, N_del, GG, Tdmax, Tdmin, Task_list_D, Cdmin, Cdmax, Rdmin, Rdmax, TaskType,  ProbabilityT, ProbabilityC, Start_tasks, Next_tasks)
        
        GG=rmnode(GG, task_for_del);
        GG_mod=GG;

        Tdmax([N_del])=[];
        Tdmin([N_del])=[];
        Cdmax([N_del])=[];
        Cdmin([N_del])=[];
        Rdmax([N_del])=[];
        Rdmin([N_del])=[];
        TaskType([N_del])=[];
        ProbabilityT([N_del])=[];
%         ProbabilityC([N_del])=[];

        Tdmax_mod=Tdmax;
        Tdmin_mod=Tdmin;
        Cdmax_mod=Cdmax;
        Cdmin_mod=Cdmin;
        Rdmax_mod=Rdmax;
        Rdmin_mod=Rdmin;
        TaskType_mod=TaskType;
        ProbabilityT_mod=ProbabilityT;
%         ProbabilityC_mod=ProbabilityC;
        Tasks_mod=string((GG_mod.Edges{:,:}));
        S_t_mod=Tasks_mod(1:end,1);
        N_t_mod=Tasks_mod(1:end,2);
        
        % ПОЧЕМУ В DEMO ВЫВОДИТ В НЕВЕРНОМ ПОРЯДКЕ/НЕ СХОДЯТСЯ С ТАБЛИЦЕЙ?
        %РЕШЕНО?
        for k=1:length(S_t_mod)
            
            ProbabilityC_mod(k)=ProbabilityC(find((S_t_mod(k)==Start_tasks)&(N_t_mod(k)==Next_tasks)));
        end

%         ismember(ismember(Start_tasks,S_t_mod),ismember(Next_tasks,N_t_mod))

    end

    function M=allcombinations(m)

        del_symb=[];
        m_buf=[];
        del_cell=[];
        [row,col]=size(m);
        for x=1:col
            del_symb=[del_symb; strip(m(1,x),"right","a")];
        end


        for x=1:col
            m_buf=[m_buf; m(:,x)];
        end
        M=[];
        M=combnk(m_buf,col);
        [row,col]=size(M);
        for x=1:col
            for y=1:row
                if count([M{y,:}],del_symb(x))>1
                    del_cell=[del_cell; y];
                end
            end
        end
        M(del_cell,:)=[];
        sort(M);
        clear m_buf; clear del_cell; clear del_symb
    end


end