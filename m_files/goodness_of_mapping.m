
Vr = [round(pk)];

B(F).noVs = RHS;
%%
timeRHSx = repmat([0, 0, 1], noLHS, 1);
timeRHSx(1:end,1:2) = RHSx;

pos = vertcat(timeRHSx,A(fr).finalg);
param.mem = 0; param.dim = 2; param.good = 2; param.quiet = 1;
res = track(pos, 3, param);

notracked = size(res,1)/2;
disp = zeros(notracked,2);

for ves = 1:2:(2*notracked);
    disp(ves,:) = res(ves,1:2)-res((ves+1),1:2);
end

dispf = disp;
dispf(~any(disp,2),:)=[];
    
sqxandy = (dispf).^2; sqdisp = sqxandy(:,1)+ sqaxandy(:,2);
avesqdisp = mean(sqdisp,1);