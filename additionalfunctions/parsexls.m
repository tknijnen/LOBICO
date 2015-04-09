function [Samples,Features,IC50s,MutationMatrix] = parsexls(file);

%% read in cls file
[num,txt,raw]=xlsread(file);
%% Parse
Samples = txt(2:end,1);
IC50s = num(:,2);

MutationMatrix = num(:,3:end);
Features = txt(1,3:end);
   
    

