function [Start_tasks, Next_tasks, ProbabilityT, ProbabilityC, Task_list_D, Tdmin, Tdmax, Cdmin, Cdmax, Rdmin, Rdmax, TaskType] = Import_table_Dh
%% Set up the Import Options and import the data
fileName = "Infocom.xlsx";
opts = spreadsheetImportOptions("NumVariables", 12);

% Specify sheet and range
opts.Sheet = "Лист2";
opts.DataRange = "A:L";

% Specify column names and types
opts.VariableNames = ["StartTasks", "NextTasks", "ProbabilityC", "TaskListD", "ProbabilityT", "Tdmin", "Tdmax", "Cdmin", "Cdmax", "Rdmin", "Rdmax", "TaskType"];
opts.VariableTypes = ["string", "string", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
% opts = setvaropts(opts, ["StartTasks", "NextTasks", "TaskListD"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["StartTasks", "NextTasks", "TaskListD"], "EmptyFieldRule", "auto");

% Import the data
tbl = readtable(fileName, opts, "UseExcel", false);
tbl(1,:)=[];

% Convert to output type
Start_tasks = tbl.StartTasks;
Next_tasks = tbl.NextTasks;
ProbabilityC = tbl.ProbabilityC;
Task_list_D = tbl.TaskListD;
% TaskType = tbl.TaskType;
% Cdmin = tbl.Cdmin;
% Cdmax = tbl.Cdmax;
% Tdmin = tbl.Tdmin;
% Tdmax = tbl.Tdmax;
% Rdmin = tbl.Rdmin;
% Rdmax = tbl.Rdmax;

Tdmin = rmmissing(tbl.Tdmin);
Tdmax = rmmissing(tbl.Tdmax);
Cdmin = rmmissing(tbl.Cdmin);
Cdmax = rmmissing(tbl.Cdmax);
Rdmin = rmmissing(tbl.Rdmin);
Rdmax = rmmissing(tbl.Rdmax);
TaskType = rmmissing(tbl.TaskType);
% Task_list_D = rmmissing(tbl.TaskListD);
ProbabilityT = rmmissing(tbl.ProbabilityT);

Task_list_D(Task_list_D=='') = [];
% TaskType(TaskType=='') = [];
% StartDate=tbl.StartDate;
% StartDate(length(Task_list_D)+1:end,:)=[];
% 
% FinishDate=tbl.FinishDate;
% FinishDate(length(Task_list_D)+1:end,:)=[];

% CT = fillmissing(tbl.CT,"constant", "FS");


end