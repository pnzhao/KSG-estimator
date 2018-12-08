function narray=binsearch2d(x,delta,varargin)
%Solve near neighbor problem using bin method. In this problem, radii are 
%different for different points.
%x:data, with the format of an [N,2] matrix.
%delta:radius 
%mx,my: number of bins in vertical and horizontal direction.
%Implementation in 2D case.
N=length(x);
b=zeros(N,2); 
%The array 'b' store the information about the index of bin that each data 
%point belongs to.
%b(1,i) and b(2,i) denotes the vertical and horizontal index of the 
%corresponding bin of x(i).
if nargin<2 || nargin>4
    error('Wrong input number');
end
if nargin==4
    mx=varargin{1};
    my=varargin{2};
else
    factor=3;
    mx=round(sqrt(N)*factor);
    my=round(sqrt(N)*factor);
end
%The above code specifies mx and my. The default value of mx and my
%increase with O(N^(-1/2)).

bincount=zeros(my,mx);
binarray=cell(my,mx);
XMIN=min(x(:,1));
XMAX=max(x(:,1));
YMIN=min(x(:,2));
YMAX=max(x(:,2));
LX=XMAX-XMIN;
LY=YMAX-YMIN;
ax=LX/mx;
ay=LY/my;
narray=zeros(1,N);
binmat=cell(my,mx);

%Calculate the number of data in each bin.
for i=1:N
    b(i,1)=ceil((x(i,1)-XMIN)/ax);
    b(i,2)=ceil((x(i,2)-YMIN)/ay);
    
    if b(i,1)==0
        b(i,1)=1;
    end
    if b(i,1)>mx
        b(i,1)=mx;
    end
    if b(i,2)==0
        b(i,2)=1;
    end
    if b(i,2)>my
        b(i,2)=my;
    end
    bincount(b(i,2),b(i,1))=bincount(b(i,2),b(i,1))+1;
    binmat{b(i,2),b(i,1)}=horzcat(binmat{b(i,2),b(i,1)},i);
end

cumcount=zeros(my,mx);
cumcount(1,1)=bincount(1,1);
for j=2:mx
    cumcount(1,j)=cumcount(1,j-1)+bincount(1,j);
end
for i=2:my
    cumcount(i,1)=cumcount(i-1,1)+bincount(i,1);
end
for i=2:my
    for j=2:mx
        cumcount(i,j)=cumcount(i-1,j)+cumcount(i,j-1)-cumcount(i-1,j-1)+bincount(i,j);
    end
end

%Begin calculating.
for i=1:N
    b0=b(i,:);
    alphax=delta(i)/ax;
    alphay=delta(i)/ay;
    
    %Now, create compare list of bins. binlist has the form [2,N],
    %and binlist(1,:) represent the x location, while binlist(2,:) 
    %represent y location. Begin with a empty list.
    binlist=[];
    
    %Consider four conditions regarding the value of alphax and alphay.
    %Note that for the first condition, there are some bins that we can
    %ensure that the distance between the query point and all the points 
    %in these bins are bounded and do not exceed delta(i), thus we do not 
    %need to calculate these distances.
    if alphax>=1 & alphay>=1
            bxplus=min(b0(1)+floor(alphax)-1,mx);
            bxminus=max(b0(1)-floor(alphax),0);
            byplus=min(b0(2)+floor(alphay)-1,my);
            byminus=max(b0(2)-floor(alphay),0);
            if bxminus==0
                c1=0;
                c2=0;
                if byminus==0
                    c3=0;
                else
                    c3=cumcount(byminus,bxplus);
                end
            else
                c2=cumcount(byplus,bxminus);
                if byminus==0
                    c1=0;
                    c3=0;
                else
                    c1=cumcount(byminus,bxminus);
                    c3=cumcount(byminus,bxplus);
                end
            end
            c4=cumcount(byplus,bxplus);

            narray(i)=c4+c1-c2-c3;
            
            %Set the boundaries.
            bx1=b0(1)-floor(alphax)-1;
            bx2=b0(1)-floor(alphax);
            bx3=b0(1)+floor(alphax);
            bx4=b0(1)+floor(alphax)+1;
            by1=max(b0(2)-floor(alphay)-1,1);
            by2=max(b0(2)-floor(alphay),1);
            by3=min(b0(2)+floor(alphay),my);
            by4=min(b0(2)+floor(alphay)+1,my);
            
            %Concatenate vertical lines.
            if bx1>=1
                binlist=horzcat(binlist,[bx1*ones(1,by4-by1+1);by1:by4]);
            end
            if bx2>=1
                binlist=horzcat(binlist,[bx2*ones(1,by4-by1+1);by1:by4]);
            end
            if bx3<=mx
                binlist=horzcat(binlist,[bx3*ones(1,by4-by1+1);by1:by4]);
            end
            if bx4<=mx
                binlist=horzcat(binlist,[bx4*ones(1,by4-by1+1);by1:by4]);
            end
            
            %Concatenate horizontal lines.
            bx5=max(bx1+2,1);
            bx6=min(bx4-2,mx);
            if by1~=b0(2)
                binlist=horzcat(binlist,[bx5:bx6;by1*ones(1,bx6-bx5+1)]);
            end
            if by4~=b0(2)
                binlist=horzcat(binlist,[bx5:bx6;by4*ones(1,bx6-bx5+1)]);
            end
            if by1~=by2
                binlist=horzcat(binlist,[bx5:bx6;by2*ones(1,bx6-bx5+1)]);
            end
            if by3~=by4
                binlist=horzcat(binlist,[bx5:bx6;by3*ones(1,bx6-bx5+1)]);
            end
    end
    
    if alphax>=1 & alphay<1
        by0=b0(2);
        by1=b0(2)-1;
        by2=b0(2)+1;
        bx1=max(b0(1)-floor(alphax)-1,1);
        bx4=min(b0(1)+floor(alphax)+1,mx);
        if by1>=1
            binlist=horzcat(binlist,[bx1:bx4;by1*ones(1,bx4-bx1+1)]);
        end
        binlist=horzcat(binlist,[bx1:bx4;by0*ones(1,bx4-bx1+1)]);
        if by2<=my
            binlist=horzcat(binlist,[bx1:bx4;by2*ones(1,bx4-bx1+1)]);
        end
    end
    
    if alphax<1 & alphay>=1
        if i==237
            a=1;
        end
        bx0=b0(1);
        bx1=b0(1)-1;
        bx2=b0(1)+1;
        by1=max(b0(2)-floor(alphay)-1,1);
        by4=min(b0(2)+floor(alphay)+1,my);
        if bx1>=1
            binlist=horzcat(binlist,[bx1*ones(1,by4-by1+1);by1:by4]);
        end
        binlist=horzcat(binlist,[bx0*ones(1,by4-by1+1);by1:by4]);
        if bx2<=mx
            binlist=horzcat(binlist,[bx2*ones(1,by4-by1+1);by1:by4]);
        end        
    end
    
    if alphax<1 & alphay<1
        bx0=b0(1);
        bx1=b0(1)-1;
        bx2=b0(1)+1;
        by0=b0(2);
        by1=max(b0(2)-1,1);
        by2=min(b0(2)+1,my);
        if bx1>=1
            binlist=horzcat(binlist,[bx1*ones(1,by2-by1+1);by1:by2]);
        end
        binlist=horzcat(binlist,[bx0*ones(1,by2-by1+1);by1:by2]);
        if bx2<=mx
            binlist=horzcat(binlist,[bx2*ones(1,by2-by1+1);by1:by2]);
        end            
    end
    
    %Now compare each bin with index stored in array binlist.
    L=length(binlist);
    for j=1:L
        B=x(binmat{binlist(2,j),binlist(1,j)},:);
        for k=1:size(B,1)
            if max(abs(x(i,:)-B(k,:)))<delta(i)
                narray(i)=narray(i)+1;
            end
        end
    end
    narray(i)=narray(i)-1; %Remove itself.
end
        