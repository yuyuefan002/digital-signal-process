%����E-M�㷨������ϸ�˹ģ��
function [im is ip times]=gmm(nim,nis,nip,data)
im=nim;
is=nis;
ip=nip;
mogs=length(ip);
l=length(data(1,:));
nocoeff=length(data(:,1));
Pwx(mogs,l)=0;
Pxw(l,mogs)=0;
dist=1;
mag=0;
dt=1;
times=0;
%ѭ��
%while dist>(mag*.01)
%while(times<time)
while and(dist>(mag*0.01),times<20)
    times=times+1;
    oldm=im;
    %E step
    %���� P(x/w)
    for i=1:l
        for j=1:mogs
            Pxw(i,j)=mvnpdf(data(:,i),im(:,j),is(:,:,j));
        end
    end
    %���� P(w/x)
    for i=1:l
	   sm=0;
         for j=1:mogs
            Pwx(j,i)=Pxw(i,j)*ip(j);
            sm=sm+Pwx(j,i);
        end
        for j=1:mogs
            Pwx(j,i)=Pwx(j,i)/sm;
        end
    end

    %M step
    %���� P(w)
    for i=1:mogs
	ip(i)=0;
        for j=1:l
            ip(i)=ip(i)+Pwx(i,j);
        end
        ip(i)=ip(i)/l;
    end

    %����Э�������
    for i=1:mogs
        sm=0;
	    is(:,:,i)=0;
        for j=1:l
            diff=data(:,j)-im(:,i);
            matr=diff*diff';
            sm=sm+Pwx(i,j);
            is(:,:,i)=is(:,:,i)+(matr.*Pwx(i,j));
        end
        is(:,:,i)=is(:,:,i)./sm;
	for k=1:nocoeff
		for m=1:nocoeff
			if k~=m
				is(k,m,i)=0;
			else
				if is(k,m,i)<0.01
					is(k,m,i)=0.01;
				end
			end
		end
	end
    end
    %�����ֵ
    for i=1:mogs
        sm=0;
	  im(:,i)=0;
        for j=1:l
            im(:,i)=im(:,i)+(data(:,j).*Pwx(i,j));
            sm=sm+Pwx(i,j);
        end
        im(:,i)=im(:,i)/sm;
    end
    dist=0;
    mag=0;
    for i=1:mogs
	tdist=0;
        tmag=0;
        for j=1:nocoeff
            tdist=tdist+(im(j,i)-oldm(j,i)).^2;
            tmag=tmag+oldm(j,i).^2;
        end
	dist=dist+sqrt(tdist);
	mag=mag+sqrt(tmag);
    end 
end
pxw=Pxw;
pwx=Pwx;
