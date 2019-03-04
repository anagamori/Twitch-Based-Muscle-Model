close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/FT';

for trialN = 300
MU_type = 'fast';

cd(data_folder)
load(['MU_' num2str(trialN)])
cd(code_folder)

[Data] = new_model_test(parameter,1,0,MU_type);
FR_half = Data{2,6};
t2t(trialN) = Data{2,5};

for i = 1:4
   if i == 1
       Lce = 0.8;
   elseif i == 2
       Lce = 0.9;
   elseif i == 3
       Lce = 1.1;
   elseif i == 4
       Lce = 1.2;
   end
   [Data] = new_model_test(parameter,Lce,FR_half,MU_type);
   
end

%close all

end