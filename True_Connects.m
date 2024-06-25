function [S_t_n, N_t_n, TrueTD, CritPath, CritTime] = True_Connects(S_t, N_t, TlD, TD)
%%
A=table(S_t, N_t);
A=unique(A,'stable');
S_t=table2array(A(:,1));
N_t=table2array(A(:,2));
S_t_n = strings([3*length(S_t),1]);
N_t_n = strings([3*length(N_t),1]);
%%
% Здесь домены равны 0
for i=1:length(S_t)
    S_t_n(i)=append(S_t(i), '_2');
    N_t_n(i)=append(N_t(i), '_1');
end
%%
% Добавляем связи для временых доменов Start_tasks
for i=length(S_t)+1:2*length(S_t)
    S_t_n(i)=append(S_t(i-length(S_t)), '_1');
    N_t_n(i)=append(S_t(i-length(S_t)), '_2');
end

% Добавляем связи для временых доменов Next_tasks
for i=2*length(S_t)+1:3*length(S_t)
    S_t_n(i)=append(N_t(i-2*length(S_t)), '_1');
    N_t_n(i)=append(N_t(i-2*length(S_t)), '_2');
end
separ=S_t_n(length(S_t)+1); %Задача определяющее начало распределения доменов времени
%%
clear A, %clear S_t, clear N_t
A=table(S_t_n, N_t_n);
A=unique(A,'stable');


S_t_n=table2array(A(:,1));
N_t_n=table2array(A(:,2));

TrueTD = zeros(length(S_t_n),1);
%Подтирание индексов для корректного сравнения
Stn=strip(S_t_n(find(S_t_n==separ):end),"right",'1');
Stn=strip(Stn,"right",'_');
%Несовпадение размерностей при данных для ганчарта
for i=find(S_t_n==separ):length(S_t_n)
    TrueTD(i)=TD(TlD==Stn(i-find(S_t_n==separ)+1));
end
clear Stn, clear A

%%
G=digraph(S_t_n, N_t_n, TrueTD);

Start_points = G.Nodes.Name(indegree(G,G.Nodes.Name)==0);
Finish_points = G.Nodes.Name(outdegree(G,G.Nodes.Name)==0);
G=addnode(G, "START");
G=addedge(G, "START", Start_points, 0);
Start_points = G.Nodes.Name(indegree(G,G.Nodes.Name)==0);

invG = G;
invG.Edges.Weight = -G.Edges.Weight;

%%


pathMatrix = cell(length(Start_points),length(Finish_points));
pathLenMatrix = zeros(length(Start_points),length(Finish_points));
for i=1:length(Start_points)
    for j =1:length(Finish_points)
%         s = Start_points{i};
%         t = Finish_points{j};
        
        [p,d] = shortestpath(invG, Start_points{i},Finish_points{j});
        if isempty(p)
            continue
        end

        pathMatrix{i,j} = p;
        pathLenMatrix(i,j) = d;
    end
end


[M,I] = min(pathLenMatrix,[],'all');
% [M,I] = max(pathLenMatrix,[],'all');
[idxrow, idxcol] = ind2sub(size(pathLenMatrix), I);

% Finish_points = pathMatrix{idxrow,idxcol};
% Finish_points=strip(Finish_points,"right",'1');
% Finish_points=strip(Finish_points,"right",'2');
% Finish_points=strip(Finish_points,"right",'_');
% Finish_points=unique(Finish_points, 'stable');
% 
% Start_points=strip(Start_points,"right",'1');
% Start_points=strip(Start_points,"right",'2');
% Start_points=strip(Start_points,"right",'_');
% Start_points=unique(Start_points, 'stable');

CritPath = pathMatrix{idxrow,idxcol};
CritPath=strip(CritPath,"right",'1');
CritPath=strip(CritPath,"right",'2');
CritPath=strip(CritPath,"right",'_');
CritPath=unique(CritPath, 'stable');
CritPath(CritPath=="START")=[];

% maxLen = pathLenMatrix(idxrow,idxcol);
CritTime = -pathLenMatrix(idxrow,idxcol);
% disp(CritTime);
end