function narray=binsearch(x,delta,m)
%Solve near neighbor problem using bin method. In this problem, radii are 
%different for different points.
%x:data.
%delta:radius 
%m:number of bins.
N=length(x);
b=zeros(1,N); 
%The array 'b' store the information about the index of bin that each data 
%point belongs to.
if nargin<2 || nargin>3
    error('Wrong input number');
end
if nargin==3
    m=varargin{1};
else
    m=round(N/10);
end
bincount=zeros(1,m);
binarray=cell(1,m,1);
MIN=min(x);
MAX=max(x);
L=MAX-MIN;
a=L/m;
narray=zeros(1,N);

%Calculate the number of data in each bin.
for i=1:N
    b(i)=ceil((x(i)-MIN)/a);
    
    %Correct abnormal values.
    if b(i)==0
        b(i)=1;
    end
    if b(i)>m
        b(i)=m;
    end
    
    bincount(b(i))=bincount(b(i))+1;
    binarray{b(i)}=[binarray{b(i)} i];
end

cumcount=zeros(1,m); %Cumulative data count.
cumcount(1)=bincount(1);
for j=2:m
    cumcount(j)=cumcount(j-1)+bincount(j);
end

%Begin calculating.
for i=1:N
    b0=b(i);
    alpha=delta(i)/a;
    
    %Determine the data that will be compared.
    if alpha>=1
        b1=floor(alpha)+b0;
        b2=ceil(alpha)+b0;
        b3=b0-floor(alpha);
        b4=b0-ceil(alpha);
        bplus=min(b0+floor(alpha)-1,m);
        bminus=max(b0-floor(alpha),0);
        if bminus==0
            narray(i)=cumcount(bplus);
        else
            narray(i)=cumcount(bplus)-cumcount(bminus);
        end
        narray(i)=narray(i)-1;
        if b1<=m
            B1=x(binarray{b1});
            for k=1:length(B1)
                if abs(x(i)-B1(k))<delta(i)
                    narray(i)=narray(i)+1;
                end
            end
        end
        if b2<=m
            B2=x(binarray{b2});
            for k=1:length(B2)
                if abs(x(i)-B2(k))<delta(i)
                    narray(i)=narray(i)+1;
                end
            end
        end
        if b3>=1
            B3=x(binarray{b3});
            for k=1:length(B3)
                if abs(x(i)-B3(k))<delta(i)
                    narray(i)=narray(i)+1;
                end
            end
        end        
        if b4>=1
            B4=x(binarray{b4});
            for k=1:length(B4)
                if abs(x(i)-B4(k))<delta(i)
                    narray(i)=narray(i)+1;
                end
            end
        end                
    else
        b1=b0+1;
        b2=b0-1;
        %Compare for B0:
        B0=x(binarray{b0});
        for k=1:length(B0)
            if abs(x(i)-B0(k))<delta(i)
                narray(i)=narray(i)+1;
            end
        end
        narray(i)=narray(i)-1; %The point compared to itself and thus add 1
        %to the result. So we need to subtract 1.
            
        %Compare for B1:
        if b1<=m
            B1=x(binarray{b1});
            for k=1:length(B1)
                if abs(x(i)-B1(k))<delta(i)
                    narray(i)=narray(i)+1;
                end
            end
        end
        
        %Compare for B2:
        if b2>=1
            B2=x(binarray{b2});
            for k=1:length(B2)
                if abs(x(i)-B2(k))<delta(i)
                    narray(i)=narray(i)+1;
                end
            end
        end
    end
end
            
            