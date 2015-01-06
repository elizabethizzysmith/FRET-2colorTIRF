Printedfilename = 'VT 40 mins 090409 analyzed 030310.tif';
TO = 76970;

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

DF = zeros(nmovies,7);
MMGR = false(512,512,nmovies);
MMVG = false(512,512,nmovies);
MM = zeros(512,512,nmovies,'uint16'); 
MM633 = zeros(512,512,nmovies,'uint16');

P(1,1:2)= [-278,-10];   %translation
P(2,1:2)= [-200,-200];   %add to redefine origin as the center of the green channel
P(3,1)= -0.003;    %rotation angle (radians)
P(4,1)= 93/100;   %scaling factor

P(7,1) = 6.4;  %Igreen bleed into red channel xtalk factor
P(7,2) = 13.5;  %Ired direct excitation of DiD factor

P(13,1) = 28; %green channel's starting point in x
P(13,2) = 477; %green channel's ending point in x

P(14,1) = 8; %green channel's starting point in y
P(14,2) = 227; %green channel's ending point in y

P(15,1) = 26; %red channel's starting point in x
P(15,2) = 508; %red channel's ending point in x

P(16,1) = 271; %red channel's starting point in y
P(16,2) = 508; %red channel's ending point in y

P(17,1) = 030210; %version of the program that was used to generate the data

P(18,1) = 103; %background in green channel w 514 
P(18,2) = 87;  %background in red channel w 514
P(19,1) = 102; %background in red channel w 633

P(20,1) = 200; %peakfind threshold for green
P(20,2) = 600; %peakfind threshold for red

P(21,1) = ((P(13,2)-P(13,1)+1))*((P(14,2)-P(14,1)+1));  %effective size of green channel in pixels
P(21,2) = (P(15,2)-P(15,1)+1)*(P(16,2)-P(16,1)+1);  %size of red channel in pixels

P(22,1) = 7566;  %514RC value that must be exceeded to be considered a fusion

for F = 1:nmovies

[A,t]=tiffreadgeneric(stackn{F});

for fr = 1:2;
    
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
  
  A(fr).imagegcf = bpass2(A(fr).imagegc,0,7);
  A(2).pksg = pkfnd(A(2).imagegcf,P(20,1),7);  
end

greens = A(2).pksg; 
notsbound = size(greens,1);
region=[-4 -4];
tallregion=repmat(region,notsbound,1);
vert=greens+tallregion;

name=strrep(stackn{F},'.tif','tsbound.rgn');
fid=fopen(name,'w');
fprintf(fid,'0 1, 1 255, 2 %u %u, 3 0 0, 4 0, 5 1, 6 2 7 7, 7 "1"\n',vert');
fclose(fid);

pk=A(1).pksr(:,1:2);
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

region=[-4 -4]; tallregion=repmat(region,2*noLHS,1); vert=Vr+tallregion;
name=strrep(stackn{F},'.tif','Vs.rgn');
fid=fopen(name,'w');
fprintf(fid,'0 1, 1 255, 2 %u %u, 3 0 0, 4 0, 5 1, 6 2 7 7, 7 "1"\n',vert');
fclose(fid);

time = datevec(A(2).datetime(12:19),13);
DF(F,1) = (3600*time(1,4) + 60*time(1,5) + time(1,6) - TO)/60; %time stamp in minutes from tiff tags and the initial time subtracted

B(F).Vs= Vr; % (x,y) positions in both green and red channels corresponding to the positions of the tethered Vs
B(F).redVs = pk; %positions on the red channel where the tethered Vs are located
B(F).greenVs = RHS;  %positions in the green channel where the tethered Vs are located
B(F).tsbound = greens; %positions in the green channel where peaks were found by searching only the green channel
B(F).notsbound = notsbound;  %number of peaks in the green channel found by searching the green channel

MM(:,:,F) = A(2).data;
MM633(:,:,F) = A(1).data;

for z = 1:notsbound
    x = greens(z,2);
    y = greens(z,1);
    MMGR(x,y,F) = 1; %tallies docking by putting ones on an matrix of zeros in the pixel corresponding the the location of the vesicles' highest intensity value
end
end

MMGR2 = MMGR; 
wr_mode = 'overwrite';

    for j = 1:nmovies
       
tmpim = MM(:,:,j);
tmpname = Printedfilename;
   
%%this loop creates a box around each V in both green and red channels,
%%which is written into the actual .tiff stack that is generated from the
%%compilation of the 514 frames (2nd frames)

peaks = double(B(j).Vs); 
pksz = size(peaks,1); 
noVs = pksz/2;

         for k = 1:pksz
            y = peaks(k,1);
            x = peaks(k,2);
            tmpim((x-3),(y-3):(y+3),1) = 65000;
            tmpim((x-2),(y+3),1) = 65000;
            tmpim((x-2),(y-3),1) = 65000;
            tmpim((x-1),(y+3),1) = 65000;
            tmpim((x-1),(y-3),1) = 65000;
            tmpim(x,(y+3),1) = 65000;
            tmpim(x,(y-3),1) = 65000;
            tmpim((x+1),(y+3),1) = 65000;
            tmpim((x+1),(y-3),1) = 65000;
            tmpim((x+2),(y+3),1) = 65000;
            tmpim((x+2),(y-3),1) = 65000;
            tmpim((x+3),(y-3):(y+3),1) = 65000;
            
            MMGR((x-3):(x+3),(y-3):(y+3),j) = 0; %blanks out the ones that correspond to the green channel's tethered V locations
         end
         
imwrite(tmpim,tmpname,'tif','Compression','none','WriteMode',wr_mode);
wr_mode = 'append';
      
B(j).noNSB = sum(sum(MMGR(:,:,j))); %finds peaks in the green channel by setting the intensities inside the boxes to zero and counting pks
 
redVs =  B(j).redVs;  
greenVs = B(j).greenVs; %positions in the red and green channel, respectively, marking locations of tethered Vs
Ired = zeros(noVs,1);  Igreen = zeros(noVs,1); Ired633 = zeros(noVs,1); Iredcorr = zeros(noVs,1); FRETeff = zeros(noVs,1);
dockscore = zeros(noVs,1); dockpos = 49*ones(noVs,1);
RedCT = P(22,1)*ones(noVs,1); %red corrected threshold with which to compared red 514 corr values to
        
 for k = 1:noVs
            y = redVs(k,1);
            x = redVs(k,2);
            yy = greenVs(k,1);
            xx = greenVs(k,2);
            
            Ired(k,1) = sum(sum(MM((x-3):(x+3),(y-3):(y+3),j)))-49*P(18,2);
            Ired633(k,1) = sum(sum(MM633((x-3):(x+3),(y-3):(y+3),j)))-49*P(19,1);
            Igreen(k,1) = sum(sum(MM((xx-3):(xx+3),(yy-3):(yy+3),j)))-49*P(18,1);
            Iredcorr(k,1) = Ired(k,1)-Igreen(k,1)/P(7,1)-Ired633(k,1)/P(7,2);
            FRETeff(k,1) = Iredcorr(k,1)/(Iredcorr(k,1)+2*Igreen(k,1));
            MMVG((xx-3):(xx+3),(yy-3):(yy+3),j)=1;
 end
 
 MMVD(:,:,j) = MMVG(:,:,j) - MMGR2(:,:,j); %to determine whether each tethered Vs has a t docked
 %to it, you subtract the array of zeros with ones at the locations
 %of the green peaks from a matrix where there are 1's in 7x7 regions
 %corresponding to the positions of the tethered Vs in the green channel
 
 for k = 1:noVs
     yy = greenVs(k,1);
     xx = greenVs(k,2);
     dockscore(k,1) = sum(sum(MMVD((xx-3):(xx+3),(yy-3):(yy+3),j)));
 end
 
 B(j).dockscore = (dockpos - dockscore);
 B(j).docked = (dockscore < dockpos); % co-localized green and red peaks
      
 %the above loop sums inside the boxes to find background-corrected
 %intensity values for red514, red633, and green514.  It calculates each
 %vesicle's "cross-talk-corrected red514 value, which is essentially the
 %intensity generated for each vesicle due to FRET. 
 
B(j).fused = (Iredcorr > RedCT);

B(j).dockorfused = B(j).fused | B(j).docked;

B(j).Ired = Ired; %list of the total intensity of each vesicle during 514 exc summed within 7x7 box (after specified bkgd subtracted)
B(j).Ired633list = Ired633;
B(j).Iredcorr = Iredcorr;
B(j).Igreen = Igreen;
B(j).FRET = FRETeff;
Int = horzcat(Igreen, Ired, Ired633, Iredcorr, FRETeff);
B(j).Int = Int;
B(j).Intdockfuse = Int(B(j).dockorfused,:);

B(j).fusedcount = sum(B(j).Intdockfuse(:,5) > 0.25);
fusedposGC = greenVs(B(j).Intdockfuse(:,5)> 0.25,:);

B(j).extra = sum(B(j).dockscore(B(j).Int(:,5) <= 0.25,:))-sum((B(j).dockorfused))+B(j).fusedcount;

B(j).tetheredVareatotal = sum(sum(MMVG(P(13,1):P(13,2),P(14,1):P(14,2),j)));

for k = 1:B(j).fusedcount
    yy = fusedposGC(k,1);
    xx = fusedposGC(k,2);
    MMVG((xx-3):(xx+3),(yy-3):(yy+3),j) = zeros;
end

B(j).tetheredVdockarea = sum(sum(MMVG(P(13,1):P(13,2),P(14,1):P(14,2),j)));
    
B(j).dockNSBcount = sum(B(j).Intdockfuse(:,5)<= 0.25) + B(j).extra;

B(j).ratioareas = B(j).tetheredVdockarea/(P(21,1)-B(j).tetheredVareatotal);
B(j).falsedock = B(j).noNSB*B(j).ratioareas;

DF(j,2) = noVs;   %total number of vesicles
%B(j).notsinUFboxes = B(j).notsbound - B(j).notsbound2;
B(j).nodockedUF = B(j).dockNSBcount - B(j).falsedock;
DF(j,3) = B(j).nodockedUF; %number of docked, unfused targets
DF(j,4) = B(j).fusedcount;   %number fused targets (which occurred after docking)
DF(j,5) = B(j).fusedcount + B(j).nodockedUF; %number docked unfused targets + number fused targets corrected for random NSB
DF(j,6) = B(j).noNSB; 
DF(j,7) = B(j).ratioareas;

% name = strrep(char(filesnames(j,1)),'.tif','Int.txt');
% dlmwrite(name,Int);

    end  
 
GDF = vertcat(B(:).Intdockfuse);
FRETeff = GDF(:,5);
name = strrep(Printedfilename,'.tif','FRETeff.txt');
dlmwrite(name, FRETeff);
    
DFcomb = sum(DF(:,2:6));
name = strrep(Printedfilename,'.tif','DFcomb.txt'); %FRET efficiency of all D&F (includes false dockings)
dlmwrite(name,DFcomb);

G = vertcat(B(:).Int); 
name=strrep(Printedfilename,'.tif','redcorr.txt'); %summary of docked (corrected for NSB) and fused (docking defined FRET <= 0.2)
dlmwrite(name,G);

name=strrep(Printedfilename,'.tif','P.txt');
dlmwrite(name,P);

name = strrep(Printedfilename,'.tif','DF.txt'); %numbers of docked and fused for each individual movie
dlmwrite(name,DF);

name = strrep(Printedfilename,'.tif','B.mat'); %saves info regarding each vesicle within each expt in a structure
save(name, 'B');
 
 









  