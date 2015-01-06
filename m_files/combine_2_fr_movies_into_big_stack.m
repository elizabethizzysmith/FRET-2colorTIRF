Printedfilename = 'D1M4 060310.tif';

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

for F = 1:nmovies
    
    [A,t]=tiffreadgeneric(stackn{F});
    MM(:,:,F) = A(2).data;
end

wr_mode = 'overwrite';

    for j = 1:nmovies

tmpim = MM(:,:,j);
tmpname = Printedfilename;
imwrite(tmpim,tmpname,'tif','Compression','none','WriteMode',wr_mode); 
  wr_mode = 'append';
    end  
    
     