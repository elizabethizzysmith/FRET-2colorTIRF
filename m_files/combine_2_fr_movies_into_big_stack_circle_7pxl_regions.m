Printedfilename = 'D1M4 071609 021510 11pxl rgns.tif';

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

MM = uint16(zeros(512,512,nmovies));
MM633 = zeros(512,512,nmovies,'uint16');

P(1,1:2)= [-275, 11];   %translation
P(2,1:2)= [-200,-200];   %add to redefine origin as the center of the green channel
P(3,1)= -0.012;    %rotation angle (radians)
P(4,1)= 91.8/100;   %scaling factor

P(7,1) = 6.4;  %Igreen bleed into red channel xtalk factor
P(7,2) = 3.9;  %Ired direct excitation of DiD factor

P(13,1) = 28; %green channel's starting point in x
P(13,2) = 477; %green channel's ending point in x

P(14,1) = 8; %green channel's starting point in y
P(14,2) = 227; %green channel's ending point in y

P(15,1) = 26; %red channel's starting point in x
P(15,2) = 508; %red channel's ending point in x

P(16,1) = 271; %red channel's starting point in y
P(16,2) = 508; %red channel's ending point in y

P(17,1) = 012710; %version of the program that was used to generate the data

P(18,1) = 105; %background in green channel w 514 
P(18,2) = 118;  %background in red channel w 514
P(19,1) = 102; %background in red channel w 633

P(20,1) = 300; %peakfind threshold for green
P(20,2) = 130; %peakfind threshold for red

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
  
end

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
Vr=[round(pk);RHS]; %vertically concatenate LHS and RHS into one big matrix - both

B(F).Vs = Vr;
B(F).redVs = pk; %positions on the red channel where the tethered Vs are located during frame 1
B(F).greenVs = RHS;  %positions in the green channel where the tethered Vs are located during frame 1

timestamp = F*ones(noLHS,1);
B(F).tethVtimest = horzcat(pk, timestamp);

MM(:,:,F) = A(2).data;
MM633(:,:,F) = A(1).data;
end
    
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
         end
         
imwrite(tmpim,tmpname,'tif','Compression','none','WriteMode',wr_mode);
wr_mode = 'append';
   
redVs =  B(j).redVs;  
greenVs = B(j).greenVs; %positions in the red and green channel, respectively, marking locations of tethered Vs
Ired = zeros(noVs,1);  Igreen = zeros(noVs,1); Ired633 = zeros(noVs,1); 

for k = 1:noVs
            y = redVs(k,1);
            x = redVs(k,2);
            yy = greenVs(k,1);
            xx = greenVs(k,2);
            
            Ired(k,1) = sum(sum(MM((x-3):(x+3),(y-3):(y+3),j)))-49*P(18,2);
            Ired633(k,1) = sum(sum(MM633((x-3):(x+3),(y-3):(y+3),j)))-49*P(19,1);
            Igreen(k,1) = sum(sum(MM((xx-3):(xx+3),(yy-3):(yy+3),j)))-49*P(18,1);
            Iredcorr(k,1) = Ired(k,1)-Igreen(k,1)/P(7,1)-Ired633(k,1)/P(7,2);
end

B(j).Iredcorr = Iredcorr;

end

name = strrep(Printedfilename,'.tif','B.mat'); %saves info regarding each vesicle within each expt in a structure
save(name, 'B');

% G = vertcat(B(:).Iredcorr); 
% name=strrep(Printedfilename,'.tif','redcorr.txt'); %summary of docked (corrected for NSB) and fused (docking defined FRET <= 0.2)
% dlmwrite(name,G);

% G = vertcat(B(:).tethVtimest); 
% name=strrep(Printedfilename,'.tif','tethVtimest.txt'); %summary of docked (corrected for NSB) and fused (docking defined FRET <= 0.2)
% dlmwrite(name,G); 



