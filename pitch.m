%������ط���ȡ��������
%clear all
function s=pitch(x,y)
%��ȡ�����źţ�ԭʼ��������Ƶ��Ϊ96000Hz
%[filename,filepath]=uigetfile('*.wav');
signal=audioread(strcat(x,y));
signal=trim(signal);
%subplot(311)
%plot(signal);
%title('ԭʼ�����ź�');
%xlabel('ʱ��')
%ylabel('����')
%signal=trim(signal);
%subplot(312)
%plot(signal);
%title('�˵���');
%xlabel('ʱ��')
%ylabel('����')


%�²�����9.6kHZ(ʵ�������źŴ�����fsͨ��ȡ7-10kHz)
%signal=decimate(signal,10);

%��֡
FrameLen = 450;
FrameInc = 128;
temp = enframe(signal, FrameLen, FrameInc);
[m,n] = size(temp);
t=1;

%ÿ֡����Ԥ����(�˲�)����ȡ��������
while t<=m
    j=1;
    for i=1:n
        newsignal(j)=temp(t,i);
        j=j+1;
    end
    
    %��ͨ�˲�
    f=900;
    fs=8000;
    [B,A]=butter(30,2*f/fs);
    newsignal=filter(B,A,newsignal);
    
    %�������޵�ƽ
    amax=max(newsignal(1:100));
    bmax=max(newsignal(n-99:n));
    cl=0.68*min(amax,bmax);
    
    %��������������ƽ����
     for i=1:n
        if newsignal(i)>cl
            newsignal(i)=newsignal(i)-cl;
        elseif newsignal(i)<-cl
            newsignal(i)=newsignal(i)+cl;
        else
            newsignal(i)=0;
        end
     end
     newsignal0=sign(newsignal);
     
     %��ȡ�źŵĻ����ֵR(k),Rk(1)��Ӧ�ڶ�ʱ����
     Rk=zeros(1,150);
     for i=21:300
         Rk(1)=Rk(1)+newsignal(i)*newsignal0(i); 
     end
     
     for i=20:150
         for j=21:300
             Rk(i)=Rk(i)+newsignal(j)*newsignal0(j+i);
         end
     end
     
     %�������ڼ�ΪʹR(k)Ϊ���ֵRmaxʱλ��k��ֵ�����Ƶ��fs�ĵ����ĳ˻�
     Rmax=max(Rk(20:150));
     if Rmax<0.25*Rk(1)
        p=0;                    %��֡Ϊ�����������������ֵΪ0
     else
        [y,p]=max(Rk(20:150));  
     end    
     s(t)=p;
     t=t+1;
end
s=s*1/fs;

%subplot(313)
%plot(s);
%title('�����켣')
%xlabel('֡��')
%ylabel('��������')