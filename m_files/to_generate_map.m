
[filename,filepath]=uigetfiles('*.tif; *.stk','Choose Movie(s)');  
filesnames=transpose(filename);
%cd(filepath); % don't need this if movies are in current dir
if iscell(filename)
   nmovies = size(filename,2);
   stackn = cell(1,nmovies);
   for i = 1:nmovies
       stackn{i} = strcat(filepath,filename{i});
   end
else
   nmovies = 1;
   stackn{1} = strcat(filepath,filename);
end

P(1,1:2)= [-278.5, -17];   %translation
P(2,1:2)= [-100,-250];   %add to redefine origin as the center of the green channel
P(3,1)= -0.003;    %rotation angle (radians)
P(4,1)= 92.1/100;   %scaling factor

P(7,1) = 6.4;  %Igreen bleed into red channel xtalk factor
P(7,2) = 15;  %Ired direct excitation of DiD factor

P(13,1) = 50; %green channel's starting point in x
P(13,2) = 466; %green channel's ending point in x

P(14,1) = 30; %green channel's starting point in y
P(14,2) = 195; %green channel's ending point in y

P(15,1) = 40; %red channel's starting point in x
P(15,2) = 490; %red channel's ending point in x

P(16,1) = 290; %red channel's starting point in y
P(16,2) = 490; %red channel's ending point in y

P(17,1) = 052510; %version of the program that was used to generate the data

P(18,1) = 100; %background in green channel w 514 
P(18,2) = 100;  %background in red channel w 514
P(19,1) = 100; %background in red channel w 633

P(20,1) = 400; %peakfind threshold for green
P(20,2) = 200; %peakfind threshold for red

for F = 1:nmovies

[A,t]=tiffreadgeneric(stackn{F});

for fr = 1;
    
  A(fr).grbkgd = P(18,1)*ones((P(13,2)-P(13,1)+1),(P(14,2)-P(14,1)+1));
  A(fr).imagegr = double(A(fr).data(P(13,1):P(13,2),P(14,1):P(14,2)))-A(fr).grbkgd;
  
  A(fr).redbkgd = P(18,2)*ones((P(15,2)-P(15,1)+1),(P(16,2)-P(16,1)+1));
  A(fr).imagered = double(A(fr).data(P(15,1):P(15,2),P(16,1):P(16,2)))-A(fr).redbkgd;  
  
  A(fr).imagegc = zeros(512);    A(fr).imagerc = zeros(512);  A(fr).image = zeros(512);    
  A(fr).imagegc(P(13,1):P(13,2),P(14,1):P(14,2)) = A(fr).imagegr;
  A(fr).imagerc(P(15,1):P(15,2),P(16,1):P(16,2)) = A(fr).imagered;
  
  A(fr).image(P(13,1):P(13,2),P(14,1):P(14,2)) = A(fr).imagegr;
  A(fr).image(P(15,1):P(15,2),P(16,1):P(16,2)) = A(fr).imagered;
   
  A(fr).imagercf = bpass2(A(fr).imagerc, 0, 7); %finds peaks in green, red channels for both frames of the movies
  A(fr).pksr = pkfnd(A(fr).imagercf,P(20,2),5); 
  A(fr).cntrpk = cntrd(double(A(fr).data), A(fr).pksr, 7);
  
  A(fr).imagegcf = bpass2(A(fr).imagegc,0,7);
  A(fr).pksg = pkfnd(A(fr).imagegcf,P(20,1),7);  
  A(fr).cntrgp = cntrd(double(A(fr).data), A(fr).pksg, 7);
end

pk=A(1).cntrpk(:,1:2);
noLHS=size(pk,1); %counts number of vesicles on the RHS (# of tethered V SNARE vesicles that were located during 633 excitation)

TL=repmat(P(1,1:2),noLHS,1);
pkRHS=pk+TL; %translation of the regions from the RHS to LHS of chip 

Odef=repmat(P(2,1:2),noLHS,1);
pkRHSO=pkRHS+Odef; %redefine the center of the LHS chip's channel as the origin

q=P(3,1);
RHSxO=P(4,1)*pkRHSO*[cos(q) -sin(q); sin(q) cos(q)]; %rotate to find real location of vesicles

UOdef=repmat(-P(2,1:2),noLHS,1);
RHSx=RHSxO+UOdef; %translate back to original coordinate system

RHS=round(RHSx); %rounds to integer values
Vr=[round(pk);RHS]; %vertically concatenate LHS and RHS into one big matrix - both channels will get regions

% Vr = round(pk);
% noVrs = size(Vr,1);

B(F).redpos = RHSx; %pos of each t from mapping the red peak to the green channel
B(F).greenpos = A(1).cntrgp(:,1:2);

region=[-4 -4]; 
tallregion=repmat(region,2*noLHS,1); 
% tallregion = repmat(region,noVrs,1);

vert=Vr+tallregion;
name=strrep(stackn{F},'.tif','map.rgn');
fid=fopen(name,'w');
fprintf(fid,'0 1, 1 255, 2 %u %u, 3 0 0, 4 0, 5 1, 6 2 7 7, 7 "1"\n',vert');
fclose(fid);

end  

redpos = vertcat(B(:).redpos); noredpos = size(redpos,1);
greenpos = vertcat(B(:).greenpos); nogreenpos = size(greenpos,1);

redmaptime = horzcat(redpos,ones(noredpos,1));
greenmaptime = horzcat(greenpos,2*ones(nogreenpos,1));

pos = vertcat(redmaptime,greenmaptime);
param.mem = 0; param.dim = 2; param.good = 2; param.quiet = 1;
res = track(pos, 3, param);

notracked = size(res,1)/2;
disp = zeros(notracked,2); veslist = zeros(notracked,2);

for ves = 1:notracked;
    vesindex = 2*ves;
    veslist(ves,:) = res(vesindex,1:2);
    disp(ves,:) = res(vesindex,1:2) - res((vesindex-1),1:2);
end

sqxandy = (disp).^2; sqdisp = sqxandy(:,1)+ sqxandy(:,2); locsqdisp = horzcat(veslist,sqdisp);
avesqdisp = mean(sqdisp,1);

z = -3:0.1:3;
subplot(2,1,1)
hist(disp(:,1),z)
title('Localization and Plotting Errors in x')
subplot(2,1,2)
hist(disp(:,2),z)
title('Localization and Plotting Errors in y')

figure2 = figure;
scatter(locsqdisp(:,1),locsqdisp(:,2),locsqdisp(:,3))
xlim(P(14,:)); ylim(P(15,:));





