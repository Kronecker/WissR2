

fid=fopen('logfile.log');

typeAlgoName={'SndRcv MasterSlave','MPI Reduce','SndRcv Tree 2 nodes', 'SndRcv Tree 16 nodes'};
procNums=[1,2,4,6,8,16,24,32];

getRMS=1;
numOfDataPoints=100;


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
        
        switch lineCell{8}
            case 'Snd/Rcv'
                switch lineCell{9}
                    case 'Tree'
                        switch lineCell{11}
                            case '2'
                                typeAlgo=3;
                            case '16'
                                typeAlgo=4;
                            otherwise
                                typeAlgo=5;
                        end
                    case 'Slave->Master'
                        typeAlgo=1;
                    otherwise
                        typeAlgo=5;
                end
            case 'MPI'
                typeAlgo=2;
            otherwise
                typeAlgo=5;
        end
        
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

errorbar(repmat(procNums',[1,numel(typeAlgoName)]),data,dataStd);

legend(typeAlgoName);
xlabel('Anzahl Prozessoren');
ylabel('Laufzeit [ms]');
