

fid=fopen('logfile.log');

typeAlgoName={'Sync'};
procNums=[1:16];

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
data=sum(data,3)./dataNums;

errorbar(repmat(procNums',[1,numel(typeAlgoName)]),data/1000,dataStd/1000);
ax=axis();
axis([ax(1:2),0,ax(4)]);
legend(typeAlgoName);
xlabel('Anzahl Prozessoren');
ylabel('Laufzeit [ms]');
