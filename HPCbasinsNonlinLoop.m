function f=HPCbasinsNonlinLoop
numsim=6;
cellOut=cell(numsim,1);
for i=1:numsim
    ci=HPCbasinsNonlin(i);
    cellOut{i}=ci;
end
save('basinsNonlinLoop.mat','cellOut')