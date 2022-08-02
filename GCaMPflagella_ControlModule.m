
% Flagella GCaMP analysis (19321)
% Select parent directory for analysis
% Input kymograps should be made by KimographClear on ImageJ
% Updated on 12/29/2021
% Updated on 5/17/2022

% Need fla_calcium6d.m



%% Parameters for kymograph analysis

% px_size = 0.104667; % pixel size in (um/px)
fps = 19.33; %frame per second
% threshold = 0.9; % between 0 and 1, larger is severe

%% Analysis and output

% select a parent directory

directory1 = dir(uigetdir);
val1 = length(directory1);
folder = directory1(1).folder;
n01 = 1;
Summary1 = table('Size',[val1, 14],'VariableTypes',{'string','int8','double','double','double','double','double','double','int8','double','double','double','double','double'});


for m01 = 1:val1
    if directory1(m01).isdir == 1 && isnumeric(regexp(directory1(m01).name, '[0-9]'))
        path1 = strcat(folder, '/', directory1(m01).name);
        disp (directory1(m01).name)
        
        cd(path1)
        list1 = ls;
        if regexp(list1, 'Summary1d') > 0
            File1 = csvread('Summary1d.csv');
            Summary1.Var1(n01) = directory1(m01).name;
            Summary1.Var2(n01) = File1(1,1);
            Summary1.Var3(n01) = File1(1,2);
            Summary1.Var4(n01) = File1(1,3);
            Summary1.Var5(n01) = File1(1,4);
            Summary1.Var6(n01) = File1(1,5);
            Summary1.Var7(n01) = File1(1,6);
            Summary1.Var8(n01) = File1(1,7);
            Summary1.Var9(n01) = File1(1,8);
            Summary1.Var10(n01) = File1(1,9);
            Summary1.Var11(n01) = File1(1,10);
            Summary1.Var12(n01) = File1(1,11);
            Summary1.Var13(n01) = File1(1,12);
            Summary1.Var14(n01) = File1(1,13);
        elseif regexp(list1, 'kymograph') > 0
            Summary1.Var1(n01) = directory1(m01).name;
            Summary2 = fla_calcium6d2(path1, fps);
            Summary1.Var2(n01) = Summary2(1,1);
            Summary1.Var3(n01) = Summary2(1,2);
            Summary1.Var4(n01) = Summary2(1,3);
            Summary1.Var5(n01) = Summary2(1,4);
            Summary1.Var6(n01) = Summary2(1,5);
            Summary1.Var7(n01) = Summary2(1,6);
            Summary1.Var8(n01) = Summary2(1,7);
            Summary1.Var9(n01) = Summary2(1,8);
            Summary1.Var10(n01) = Summary2(1,9);
            Summary1.Var11(n01) = Summary2(1,10);
            Summary1.Var12(n01) = Summary2(1,11);
            Summary1.Var13(n01) = Summary2(1,12);
            Summary1.Var14(n01) = Summary2(1,13);
        end
        n01 = n01+1;
        cd(folder)
    end
    
end

% Remove empty rows
x1 = ismissing(Summary1.Var1);
Summary1(x1,:) = [];

cd(folder)
writetable(Summary1,'Summary_all1d.txt','Delimiter','\t');