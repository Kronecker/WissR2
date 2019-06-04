

fid=fopen('logfile.log');

typeAlgoName={'Sync'};
procNums=[1:16,32,48,64];

getRMS=1;
numOfDataPoints=10;


data=zeros([numel(procNums),numel(typeAlgoName),numOfDataPoints]);
dataNums=zeros([numel(procNums),numel(typeAlgoName)]);
dataStd=zeros([numel(procNums),numel(typeAlgoName)]);



try
    while(~feof(fid))
        line=fgetl(fid);
        
        lineCell=strsplit(line,' ');
        procs=str2double(lineCell{2});
        procId=find(procNums==procs,1);
        if(isempty(procId))
           continue; 
        end
        
        times=str2double(lineCell{6}(1:end-2));
        
        typeAlgo=1;
        
        if(getRMS)
            dataNums(procId,typeAlgo)=dataNums(procId,typeAlgo)+1;
            data(procId,typeAlgo,dataNums(procId,typeAlgo))=times;
        else
            data(procId,typeAlgo)=data(procId,typeAlgo)+times;
            dataNums(procId,typeAlgo)=dataNums(procId,typeAlgo)+1;
        end
        
    end
    
    
catch me
    fclose(fid);
    rethrow(me);
end
fclose(fid);

dataStd=std(data,0,3);
dataSum=sum(data,3)./dataNums;


data=squeeze(data);
data=sort(data,2);

data=data(:,1:end-3);
dataNumsCut=dataNums-3;

dataStd=std(data,0,2);

dataSum=sum(data,2)./dataNumsCut;

errorbar(repmat(procNums(1:end-3)',[1,numel(typeAlgoName)]),(dataSum(1:end-3))/1000,dataStd(1:end-3)/1000);
ax=axis();
axis([0.75 ax(2),0,ax(4)]);
legend(typeAlgoName);
xlabel('Anzahl MPI Nodes');
ylabel('Laufzeit [ms]');
title([{'2000 Jacobi-Iteration auf 1024x1024 inneren Gitterpunkten'},{'Asynchrone Kommunikation mit MPI'},{'CPU: 2 physische Kerne mit je 2 logischen Kernen'}])
