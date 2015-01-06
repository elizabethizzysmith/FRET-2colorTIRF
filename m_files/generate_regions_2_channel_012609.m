%Choose the movies to analyze%% Define the parameters 

P(1,1:2)= [-274, -7];   %translation
P(2,1:2)= [-100,-250];   %add to redefine origin as the center of the green channel
P(3,1)= -0.003;    %rotation angle (radians)
P(4,1)= 92.1/100;   %scaling factor

P(7,1) = 6.4;  %Igreen bleed into red channel xtalk factor
P(7,2) = 15;  %Ired direct excitation of DiD factor

P(13,1) = 50; %green channel's starting point in x
P(13,2) = 466; %green channel's ending point in x

P(14,1) = 23; %green channel's starting point in y
P(14,2) = 207; %green channel's ending point in y

P(15,1) = 40; %red channel's starting point in x
P(15,2) = 490; %red channel's ending point in x

P(16,1) = 290; %red channel's starting point in y
P(16,2) = 490; %red channel's ending point in y

P(17,1) = 052510; %version of the program that was used to generate the data

P(18,1) = 100; %background in green channel w 514 
P(18,2) = 100;  %background in red channel w 514
P(19,1) = 100; %background in red channel w 633

P(20,1) = 700; %peakfind threshold for green
P(20,2) = 200; %peakfind threshold for red


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

F=1:nmovies; 
B(F).filename=stackn{F};

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
   
  A(fr).imagercf = bpass2(A(fr).imagerc, 0, 5); %finds peaks in green, red channels for both frames of the movies
  A(fr).pksr = pkfnd(A(fr).imagercf,P(20,2),5); 
  
  A(fr).imagegcf = bpass2(A(fr).imagegc,0,5);
  A(fr).pksg = pkfnd(A(fr).imagegcf,P(20,1),5);  
  
end

pk = A(fr).pksg;

nopks = size(pk,1);
Vr33=repmat(pk,1,33);
region=[-5	-5	-3	4	-2	4	-1	4	0	4	0	3	1	3	2	2	2	1	3	1	3	0	3	-1	3	-2	2	-2	2	-3	1	-4	0	-4	0	-5	-1	-5	-2	-5	-3	-5	-3	-4	-4	-4	-5	-3	-5	-2	-6	-2	-6	-1	-6	0	-6	1	-5	1	-5	2	-4	3	-3	3];
tallregion=repmat(region,nopks,1);
vert=Vr33+tallregion;

name=strrep(stackn{F},'.tif','greenpks.rgn');
fid=fopen(name,'w');
fprintf(fid,'0 3, 1 255, 2 %u %u, 3 0 0, 4 0, 5 1, 6 32 %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u, 7 "1"\n',vert')
fclose(fid);







