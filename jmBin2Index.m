function index = jmBin2Index(bin,thr)

n=1;index=[];

    for i=2:size(bin,1)
        if bin(i-1)<thr && bin(i)>=thr
            index(n,1)=i;
        elseif bin(i-1)>=Params.NumTTthreshold && bin(i)<thr
            index(n,2)=i;
            n=n+1;
        end
    end
    
    if bin(1)>=thr
        index(1,1)=1;
    end