
% 5/14/2021
% 9/6/2021 Updated
% Control script for detailed analysis of calcium influx and IFT
% Output 
% 12/14/2021 Updated
% 12/17/2021 Updated
% 2/3/2022 Updated
% 5/27/2022 Updated

% Need following finctions: fla_calciumD5b.m, formatKymoButlerC2.m, findPeaks2.m

%% Parameters

fps = 9.66; % frames per second
px_size = 0.104667; % pixel size in (um/px)


%% 

directory01 = dir(uigetdir);

val1 = size(directory01,1);
folder = directory01(1).folder;
n01 = 1;
Summary1 = table('Size',[val1, 47],'VariableTypes',{'string','uint16','double','double','double','double','double','double','uint16','double','double','double','double','double','double','double','double','double','double','uint16','double','double','double','double','double','double','uint16','double','double','double','double','double','double','uint16','double','double','double','double','double','double','uint16','double','double','double','double','double','double'});
% Need to fix variabletypes!!!


for m01 = 1:val1
    close
    if directory01(m01).isdir == 1 && ~isempty(regexp(directory01(m01).name, '[0-9]', 'once'))
        path1 = strcat(folder, '/', directory01(m01).name);
        disp (directory01(m01).name)
           
        cd(path1)
        list1 = ls;
        
    [IFT2, Cal2] = fla_calciumD5f(path1, px_size, fps);
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
        Summary1.Var13(n01) = IFT2(1,12);
        Summary1.Var14(n01) = IFT2(1,13);
        Summary1.Var15(n01) = IFT2(1,14);
        Summary1.Var16(n01) = Cal2(1,1);
        Summary1.Var17(n01) = Cal2(1,2);
        Summary1.Var18(n01) = Cal2(1,3);
        Summary1.Var19(n01) = Cal2(1,4);
        Summary1.Var20(n01) = Cal2(1,5);
        Summary1.Var21(n01) = Cal2(1,6);
        Summary1.Var22(n01) = Cal2(1,7);
        Summary1.Var23(n01) = Cal2(1,8);
        Summary1.Var24(n01) = Cal2(1,9);
        Summary1.Var25(n01) = Cal2(1,10);
        Summary1.Var26(n01) = Cal2(1,11);
        Summary1.Var27(n01) = Cal2(1,12);
        Summary1.Var28(n01) = Cal2(1,13);
        Summary1.Var29(n01) = Cal2(1,14);
        Summary1.Var30(n01) = Cal2(1,15);
        Summary1.Var31(n01) = Cal2(1,16);
        Summary1.Var32(n01) = Cal2(1,17);
        Summary1.Var33(n01) = Cal2(1,18);
        Summary1.Var34(n01) = Cal2(1,19);
        Summary1.Var35(n01) = Cal2(1,20);
        Summary1.Var36(n01) = Cal2(1,21);
        Summary1.Var37(n01) = Cal2(1,22);
        Summary1.Var38(n01) = Cal2(1,23);
        Summary1.Var39(n01) = Cal2(1,24);
        Summary1.Var40(n01) = Cal2(1,25);
        Summary1.Var41(n01) = Cal2(1,26);
        Summary1.Var42(n01) = Cal2(1,27);
        Summary1.Var43(n01) = Cal2(1,28);
        Summary1.Var44(n01) = Cal2(1,29);
        Summary1.Var45(n01) = Cal2(1,30);
        Summary1.Var46(n01) = Cal2(1,31);
        Summary1.Var47(n01) = Cal2(1,32);
    end
n01 = n01+1;
cd(folder)
end
            
% Remove empty rows
x1 = ismissing(Summary1.Var1);
Summary1(x1,:) = [];

% Export summary data
cd(folder)
writetable(Summary1,'Summary_allD5f.txt','Delimiter','\t');

