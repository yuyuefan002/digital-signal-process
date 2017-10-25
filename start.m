
close all;
chos=0;
possibility=8;

messaggio='Insert the number of set: each set determins a class. This set should include a number of speech for each person, with some variations in expression and in the lighting.';

while chos~=possibility,
    chos=menu('������Ƶ�����봦��ϵͳ','ѡ����Ƶ����������ݿ�','ѡ����Ƶ����������ʶ��','ɾ�����ݿ��е�����',...
        'ѡ����Ƶ���з���','˫������Ƶ��Ϊ������','����ȥ��','��Դ����','�˳�');
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % ��һѡ��洢�����Ƶ����
    if chos==1
        clc;
        close all;
        selezionato=0;
        while selezionato==0                                                %ֱ��ѡ���ļ�ֹͣ
            [namefile,pathname]=uigetfile({'*.wav','speech Files (*.wav)'},' ѡ����Ƶ�ļ�');
            if namefile~=0
     
                selezionato=1;
            else
                disp('ѡ����Ƶ�ļ�');
            end
           
        end
        filt=melfilter(150,300,15);                                         %���Ǵ�ͨ�˲���1.��Ƶ�׽���ƽ�� 2.������Ϣ�������ƥ���ٶ�
        fr1=frm(strcat(pathname,namefile),16,8000,1);                       %���룬��֡
        mc2=train(fr1,filt,20);                                             %��ȡ÷��Ƶ�ʵ���ϵ��������ٲ�ͼ
        mc2=mc2(3:18,:);
        mc1=banshengsin(mc2);                                               %Ϊ����ϵ���˰������Һ���
        s1=pitch(pathname,namefile);                                        %��ȡ��������
        a=length(s1);
        b=length(mc1(1,:));
        if a>b                                                              
            s1(b+1:a)=[];                                                   %�������ڳ�������Ĳ�����Ϊ��
        else
            s1(a+1:b)=0;                                                    %�������ڶ̣���ô����һ����Ϊ0
        end
        mc1=[mc1;s1];
        [im is ip]=init(mc1,16);                                            %������ϸ�˹ģ�͵ĳ�ʼֵ����ֵ�����ϵ��
        [nim1 nis1 nip1 times]=gmm(im,is,ip,mc1);                           %����E-M�㷨������ϸ�˹ģ��
      data=struct('name',{},'means',{},'cov',{},'prob',{},'pitch',{});
        
        if (exist('speech_database.dat')==2)
           load('speech_database.dat','-mat');
            speaker_number=speaker_number+1;
           prompt={'��������������'};
   name='�������� ';
   numlines=1;
   defaultanswer={'��'};
   answer=inputdlg(prompt,name,numlines,defaultanswer);
   data(speaker_number).name=answer{1,1};                                   %�����û���Ϣ
          data(speaker_number).means=nim1;
          data(speaker_number).cov=readcov(nis1);
          data(speaker_number).prob=nip1;
          data(speaker_number).pitch=s1;
          save('speech_database.dat','data','speaker_number','-append');
       else
           speaker_number=1;   
    prompt={'�������û�������'};
   name='�û� ';
   numlines=1;
   defaultanswer={'��'};
   answer=inputdlg(prompt,name,numlines,defaultanswer);                     %������������
   data(speaker_number).name=answer{1,1};                                   %����������ִ���data���ݱ���
           data(speaker_number).means=nim1;
        data(speaker_number).cov=readcov(nis1);
        data(speaker_number).prob=nip1;
        data(speaker_number).pitch=s1;
           save('speech_database.dat','data','speaker_number');
       end
        
        message=strcat(answer{1,1},'����Ϣ�ѳɹ��������ݿ��� ')
        msgbox(message,'speechsignal DataBase','help')
    end
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %�ڶ�ѡ�����ƥ��
    if chos==2
        clc;
        close all;
        selezionato=0;
        while selezionato==0
            [namefile,pathname]=uigetfile({'*.wav','speech Files (*.wav)'},'Chose speech signal');
            if namefile~=0
             
                selezionato=1;
            else
                disp('Select a speech signal');
            end
          
        end
      
        if (exist('speech_database.dat')==2)
            load('speech_database.dat','-mat');
            filt=melfilter(150,300,15);
 fr=frm(strcat(pathname,namefile),16,8000,3);                               %��ȡ��Ƶ���н���
 l=length(fr(1,:));                                                         %��ȡ��һ�еĳ���
 nosp=length(data);
 k=0;
 b=0;
 r=nosp;
 while(r~=1)                                                                %��k����data�е���������2�Ķ��ٴη�
  r=floor(r/2);
  k=k+1;
end
p(2,nosp)=0;p(1,1)=0;
for i=1:nosp
    p(2,i)=i;
end
mc4=train(fr,filt,20);                                                      %��ȡ÷��Ƶ�ʵ���ϵ��������ٲ�ͼ
mc4=mc4(3:18,:);
mc=banshengsin(mc4);                                                        %Ϊ����ϵ���˰������Һ���
pitch2=pitch(pathname,namefile);                                            %��ȡ��������
a=length(pitch2);
b=length(mc(1,:));
if a>b                                                                      %�Ƚϻ������ں�˵�������ĳ��ȣ��������ڳ�����ôȥ���������ڶ�����ֵ�����
    pitch2(b+1:a)=[];
else
    pitch2(a+1:b)=0;                                                        %���˵���������ȳ���Ϊ��������ȱ�ٵĳ��Ȳ�0
end
mc=[mc;pitch2];                                                             %���������ڲ������һ��
coff=length(mc(:,1));                                                       %��ȡ��һ�г���
o=length(mc(1,:));                                                          %��ȡ��һ�еĳ���
frameparts=struct('frame',{});
s=mod(l,k);
y=floor(l/k);
if s==0
   for i=1:k                                                                %data�ĳ�����2��k�η�
       frameparts(i).frame(coff,y)=0;                                       %���������Ƶ������Ϊ������yΪ�����������һ��ֵ��0
   end                                                                      %��Ҫ����Ϊ������frameparts�ľ����
else
    for i=1:s
      frameparts(i).frame(coff,y+1)=0;
    end
    for i=s+1:k
      frameparts(i).frame(coff,y)=0;
    end
end
for r=1:k
 count=1;
   for i=r:k:l
      frameparts(r).frame(:,count)=mc(:,i);
      count=count+1;
   end
end
c=length(data);                                                             %���ݱ����г�Ա�ĸ���
for  i=1:k
   % tic
   p1=ident2(frameparts(i).frame,filt,data,p);                              %���������Ȼ����Ѱ���������ݵĸ�˹ģ��
 %  toc
   p=upd_pr(p,p1);                                                          %��������p1ȥ����p
   p=nmax1(p);                                                              %ȡp�����ֵ
end
p2=p(1)/o;
scores=zeros(nosp,1);
for i=1:nosp                                                                %nosp�������ݵĳ���
   pitch1=data(i).pitch';
  % tic
   scores(i,1)=myDTW(pitch2,pitch1(1:length(pitch2)));                      %���ر����������������ݵ����·��
  % toc
end
scores;
[m,n]=sort(scores);                                                         %��scores��������

b=p(2,1);
if or((p2>-25),b==n)
nm=data(b).name;
       message=strcat('˵�����ǣ� ',nm);
       msgbox(message,'DataBase Info','help');
else
    message='�޷�ƥ�䵽�����Ա';
    msgbox(message,'DataBase Info','help');
end
        else
            message='���ݱ�Ϊ�գ�������û�';
            msgbox(message,'speech DataBase Error','warn');    
        end
        
    end 
   
%����ѡ�ɾ��ȫ�����ݣ���ֻɾ��һ���˵�����    
    if chos==3                                                  
        clc;
        close all;
        if (exist('speech_database.dat')==2)                                %�����б�����κ��˵������ͻ᷵��2
             load('speech_database.dat','-mat');                            %�����ݱ���matlab��׼��ʽ��
            button = questdlg('��������ɾ����һλ������?',...
                   'ѡ��',...
                   '������','ָ����һλ','all');
             if strcmp(button,'������')                                     %�жϰ����ĸ���ť
                delete('speech_database.dat');
                msgbox('���ѳɹ�ɾ����������','Database removed','help');
             else 
                  prompt={'����������ɾ�����ݵ���Ա����'};
   name='ɾ��ָ����Ա';
   numlines=1;
   defaultanswer={'��'};
   answer=inputdlg(prompt,name,numlines,defaultanswer);                     %���������answerΪ���������
nspeaker=length(data)                   ;
 names=cell(1,nspeaker);
 for i=1:nspeaker                                                           %����һ�����ڼ���������
names{1,i}=data(i).name;
 end
[a,b]=ismember(answer{1,1},names);                                          %�ж�����������Ƿ��ڱ���                            
 if a==0
       warndlg('������ָ����ѡ','���棡')
   else 
      data(b)=[];                                                           %�ҵ��˸�������Ա����λ�ã���ոô�����Ƶ
      speaker_number=length(data);
      save('speech_database.dat','data','speaker_number','-append');
   message=strcat('�ɹ�ɾ�� ',answer{1,1},'����Ϣ');
            msgbox(message,'�ɹ�ɾ��','help');
 end
             end
        else
            warndlg('���ݱ���Ϊ��.',' ���棡 ')
        end
    end 
    %����ѡ���Ƶ�źŷ���
    if chos==4
        clc;
        close all;
        selezionato=0;
        while selezionato==0                                                %û��ѡ���ļ���һֱѡ��
            [namefile,pathname]=uigetfile({'*.wav','speech signal (*.wav)'},'��ѡ����Ƶ�ļ�');
            if namefile~=0
               [x,fs]=audioread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('��ѡ����Ƶ�ļ�');
            end
        
        end
        %����ͼ���ٲ�ͼ��
        figure('Name','�ٲ�ͼ');
        Y1 = x(:,1);        %YΪ˫�������ݣ�ȡ��2ͨ��
        plot(Y1);           %��Y1����ͼ
        grid on;
        specgram(Y1,2048,44100,2048,1536);
         %Y1Ϊ��������
         %FFT֡��2048��(��44100HzƵ��ʱԼΪ46ms)
         %����Ƶ��44.1KHz
         %�Ӵ����ȣ�һ����֡�����
         %֡�ص����ȣ��˴�ȡΪ֡����3/4
         T=1/fs;
         N=length(x);
         n=0:N-1;
         t=n/fs;
         figure('Name','Ƶ��');
         subplot(2,2,1);plot(t,x);axis([0 N/fs -0.4 0.4]);
         title('ʱ���ź�');
         ylabel('���');xlabel('t');
         N0=256;
         n0=0:N0-1;
         f=n*fs/N;
         z=fft(x,N);mf1=abs(z');
         subplot(2,2,2);plot(f,mf1);xlabel('f');ylabel('��ֵ');
         subplot(2,2,3);plot(f,angle(z));xlabel('��λ');ylabel('��ֵ');
         subplot(2,2,4);plot(freqz(x));


        
   %[im is ip]=init(mc1,16);
    % [nim1 nis1 nip1 times]=gmm(im,is,ip,mc1);
     
    end  
         if chos==5
        clc;
        close all;
        selezionato=0;
        while selezionato==0                                                %û��ѡ���ļ���һֱѡ��
            [namefile,pathname]=uigetfile({'*.wav','speech signal (*.wav)'},'��ѡ����Ƶ�ļ�');
            if namefile~=0
               [x,fs]=audioread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('��ѡ����Ƶ�ļ�');
            end
        
        end
            x=x(:,1);
            audiowrite(namefile,x,fs);
            
         end
         %�������֣�����ȥ��
              if chos==6
        clc;
        close all;
        selezionato=0;
        while selezionato==0                                                %û��ѡ���ļ���һֱѡ��
            [namefile,pathname]=uigetfile({'*.wav','speech signal (*.wav)'},'��ѡ����Ƶ�ļ�');
            if namefile~=0
               [sound,fs]=audioread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('��ѡ����Ƶ�ļ�');
            end
        
        end
            sound=sound(:,1);                                               %���뵥������Ƶ�ź�
            count=length(sound);
            t=(0:count-1)/fs;                                               %�ֶδ���
            sound=sound';
            [c,l]=wavedec(sound,3,'db6');                                   %��С������db6���źŽ�������ֽ�
            sigma=wnoisest(c,l,1);                                          %���Ƴ߶�Ϊ1��������׼ƫ��
            alpha=2;                                                        %����������һ��ȡ2
            thr=wbmpen(c,l,sigma,alpha);                                    %��ȡ��������е���ֵ
            keepapp=1;  
            yd=wdencmp('gbl',c,l,'db6',3,thr,'s',keepapp);                  %���źŽ�������
            subplot(1,2,1)
            plot(t,sound);
            title('ԭʼ�����ź�');
            axis([0 11 min(sound)-0.1 max(sound)+0.1])
            subplot(1,2,2)
            plot(t,yd);
            title('ȥ�������ź�');
            axis([0 11 min(yd)-0.1 max(yd)+0.1])
            audiowrite('final.wav',yd,fs);
              end
              %���߲��֣���Դ����
              if chos==7
        clc;
        close all;
        selezionato=0;
        while selezionato==0                                                %û��ѡ���ļ���һֱѡ��,
            [namefile,pathname]=uigetfile({'*.wav','speech signal (*.wav)'},'��ѡ���һ����Ƶ�ļ�');
            if namefile~=0
               [I1,fs]=audioread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('��ѡ����Ƶ�ļ�');
            end
        
        end
        selezionato=0;  
        while selezionato==0                                                %û��ѡ���ļ���һֱѡ��
            [namefile,pathname]=uigetfile({'*.wav','speech signal (*.wav)'},'��ѡ��ڶ�����Ƶ�ļ�');
            if namefile~=0
               [I2,fs]=audioread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('��ѡ����Ƶ�ļ�');
            end
        
        end
        selezionato=0;  
        while selezionato==0                                                %û��ѡ���ļ���һֱѡ��
            [namefile,pathname]=uigetfile({'*.wav','speech signal (*.wav)'},'��ѡ����Ķ���Ƶ�ļ�');
            if namefile~=0
               [I3,fs]=audioread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('��ѡ����Ƶ�ļ�');
            end
        
        end
        subplot(4,3,1),plot(I1),title('�����ź�1');
        subplot(4,3,2),plot(I2),title('�����ź�2');
        subplot(4,3,3),plot(I3),title('�����ź�3');


        l1=length(I1);
        l2=length(I2);
        l3=length(I3);
        
        if l1>=l2&&l1>=l3                                                %ά��ͳһ��
             I2(l2:l1)=0;
             I3(l3:l1)=0;
        end
        if l2>=l1&&l2>=l3
              I1(l1:l2)=0;
              I3(l2:l2)=0;
        end
        if l3>=l1&&l3>=l2
              I1(l1:l3)=0;
              I2(l2:l3)=0;
        end
        II1=I1';
        II2=I2';
        II3=I3';
        S=[II1;II2;II3];
        Sweight=rand(size(S,1));
        MixedS=Sweight*S;
                                                                         %����Ͼ����������в����
        subplot(4,3,4),plot(MixedS(1,:)),title('����ź�1');
        subplot(4,3,5),plot(MixedS(2,:)),title('����ź�2');
        subplot(4,3,6),plot(MixedS(3,:)),title('����ź�3');
       
        audiowrite('1mixwav1.wav',MixedS(1,:),fs);                       %����wav����
        audiowrite('1mixwav2.wav',MixedS(2,:),fs);
        audiowrite('1mixwav3.wav',MixedS(3,:),fs);
        MixedS_bak=MixedS;
%%%%%%%%%%%%%%%%%%%%%%%%%%��׼��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ȥ��ֵ��ʹ����Ϊ���ֵ
        MixedS_mean=zeros(3,1);
        for i=1:3                                                        %����MixedS�ľ�ֵ
            MixedS_mean(i)=mean(MixedS(i,:));
        end                                                      

        for i=1:3                                                        %�Ծ�����ÿһ��ֵ��ȥ��ֵ��ʹ�þ����Ϊ���ֵ
             for j=1:size(MixedS,2)
                MixedS(i,j)=MixedS(i,j)-MixedS_mean(i);
             end
        end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%�׻�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ȥ������ԣ����Խ�������ĸ����ԣ�ʹ���㷨�����Ա�ã��򻯺���������������ȡ����
            
        MixedS_cov=cov(MixedS');                                         %covΪ��Э����ĺ���
        [E,D]=eig(MixedS_cov);                                           %���źž����Э�������������ֵ�ֽ�
        %Q=inv(sqrt(D))*(E)';                                            %QΪ�׻�����
        Q=sqrt(D)\(E)';
        MixedS_white=Q*MixedS;                                           %MixedS_whiteΪ�׻�����źž���
        IsI=cov(MixedS_white');                                          %IsIӦΪ��λ��
	
%%%%%%%%%%%%%%%%%%%%%%%%FASTICA�㷨%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%������Ѱ����Դ�ķǸ�˹���Ĳ������������
        X=MixedS_white;                                                  %�����㷨����X���в���
        [VariableNum,SampleNum]=size(X);                                 
        numofIC=VariableNum;                                             %�ڴ�Ӧ���У�����Ԫ�������ڱ�������(����Դ����)
        B=zeros(numofIC,VariableNum);                                    %��ʼ��������w�ļĴ����,B=[b1 b2 ... bd]
        for	r=1:numofIC                                                  %������Դ�ĸ��������Ӧ�����Ķ�����Դ
            i=1;maxIterationsNum=100;                                    %����������������������ÿ�������������Ե������������˴�����
            IterationsNum=0;	
            b=rand(numofIC,1)-.5;                                        %�������b��ֵ,
            b=b/norm(b);                                                 %��b��׼��norm(b):����Ԫ��ƽ���Ϳ�����,��b���й�һ��
            while i<=maxIterationsNum+1                                  %���ж�ε������ƽ��������������
                if i == maxIterationsNum                                 %�����������ڲ�������ѭ����������
                        fprintf('\n��%d������%d�ε����ڲ���������', r,maxIterationsNum); 
                        break	;
                end	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���Ĺ�ʽ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                bOld=b;                                 
                a2=1;	
                u=1;	
                t=X'*b;	
                g=t.*exp(-a2*t.^2/2);	
                dg=(1-a2*t.^2).*exp(-a2*t.^2/2);	
                b=((1-u)*t'*g*b+u*X*g)/SampleNum-mean(dg)*b;	
                b=b-B*B'*b;                                              %��b������
                b=b/norm(b);                                             %��b��һ��
                if abs(abs(b'*bOld)-1)<1e-9                              %����㹻���ƣ���
                        B(:,r)=b;                                        %������������b
                        break;
                end 
            i=i+1; 
             end	
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%ICA��������ݸ�ԭ����ͼ%%%%%%%%%%%%%%%%%%%%%%%%%
    ICAedS=B'*Q*MixedS_bak;                                         %����ICA��ľ��󣬵õ����ź������Դ�źŽ��зֽ�
                                                                    %����Ͼ����������в����
    subplot(4,3,7),plot(ICAedS(1,:)),title('ICA����ź�1');
    subplot(4,3,8),plot(ICAedS(2,:)),title('ICA����ź�2');
    subplot(4,3,9),plot(ICAedS(3,:)),title('ICA����ź�3');

    audiowrite('1dstwav1.wav',ICAedS(1,:),fs);                      %����wav����
    audiowrite('1dstwav2.wav',ICAedS(2,:),fs);
    audiowrite('1dstwav3.wav',ICAedS(3,:),fs);

              end
end