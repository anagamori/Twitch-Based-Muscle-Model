
%rng(1)
N_MU = 200;
CT = raylrnd(20,1,N_MU)+20;
CT = sort(CT,'descend');
figure(2)
histogram(CT)
min(CT)
max(CT)

%%
cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
save('CT_target','CT')
cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code')