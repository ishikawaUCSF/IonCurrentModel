
% 8/15/2021
% 5/2/2022 Updated
% 5/4/2022 Updated
% Control script for detailed analysis of KAP-GFP kymographs


% Need following finctions: "KinesinKymoButler2.m", formatKymoButlerC1.m, 

%% Parameters

fps = 19.33; % frames per second
px_size = 0.104667; % pixel size in (um/px)


%% Analysis

directory01 = dir(uigetdir);

val1 = size(directory01,1);
folder = directory01(1).folder;
n01 = 1;
Summary1 = table('Size',[val1, 12],'VariableTypes',{'string','uint16','double','double','double','double','double','double','double','double','double','double'});

for m01 = 1:val1
    close
    if directory01(m01).isdir == 1 && ~isempty(regexp(directory01(m01).name, '[0-9]', 'once'))
        path1 = strcat(folder, '/', directory01(m01).name);
        disp (directory01(m01).name)
           
        cd(path1)
        list1 = ls;
        
    IFT2 = KinesinKymoButler2(path1, px_size, fps);
        Summary1.Var1(n01) = directory01(m01).name;
        Summary1.Var2(n01) = IFT2(1,1);
        Summary1.Var3(n01) = IFT2(1,2);
        Summary1.Var4(n01) = IFT2(1,3);
        Summary1.Var5(n01) = IFT2(1,4);
        Summary1.Var6(n01) = IFT2(1,5);
        Summary1.Var7(n01) = IFT2(1,6);
        Summary1.Var8(n01) = IFT2(1,7);
        Summary1.Var9(n01) = IFT2(1,8);
        Summary1.Var10(n01) = IFT2(1,9);
        Summary1.Var11(n01) = IFT2(1,10);
        Summary1.Var12(n01) = IFT2(1,11);
    end
n01 = n01+1;
cd(folder)
end
            
% Remove empty rows
x1 = ismissing(Summary1.Var1);
Summary1(x1,:) = [];

% Export summary data
cd(folder)
writetable(Summary1,'Summary_allC2.txt','Delimiter','\t');

