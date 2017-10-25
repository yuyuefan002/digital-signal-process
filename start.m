
close all;
chos=0;
possibility=8;

messaggio='Insert the number of set: each set determins a class. This set should include a number of speech for each person, with some variations in expression and in the lighting.';

while chos~=possibility,
    chos=menu('数字音频分析与处理系统','选择音频并添加至数据库','选择音频并进行声纹识别','删除数据库中的内容',...
        '选择音频进行分析','双声道音频改为单声道','语音去噪','声源分离','退出');
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % 第一选项：存储相关音频特征
    if chos==1
        clc;
        close all;
        selezionato=0;
        while selezionato==0                                                %直到选中文件停止
            [namefile,pathname]=uigetfile({'*.wav','speech Files (*.wav)'},' 选择音频文件');
            if namefile~=0
     
                selezionato=1;
            else
                disp('选择音频文件');
            end
           
        end
        filt=melfilter(150,300,15);                                         %三角带通滤波，1.对频谱进行平滑 2.降低信息量，提高匹配速度
        fr1=frm(strcat(pathname,namefile),16,8000,1);                       %降噪，分帧
        mc2=train(fr1,filt,20);                                             %提取梅尔频率倒谱系数，组成瀑布图
        mc2=mc2(3:18,:);
        mc1=banshengsin(mc2);                                               %为输入系数乘半升正弦函数
        s1=pitch(pathname,namefile);                                        %获取基因周期
        a=length(s1);
        b=length(mc1(1,:));
        if a>b                                                              
            s1(b+1:a)=[];                                                   %基音周期长，多出的部分设为空
        else
            s1(a+1:b)=0;                                                    %基音周期短，那么增加一部分为0
        end
        mc1=[mc1;s1];
        [im is ip]=init(mc1,16);                                            %给出混合高斯模型的初始值：均值，方差，系数
        [nim1 nis1 nip1 times]=gmm(im,is,ip,mc1);                           %利用E-M算法建立混合高斯模型
      data=struct('name',{},'means',{},'cov',{},'prob',{},'pitch',{});
        
        if (exist('speech_database.dat')==2)
           load('speech_database.dat','-mat');
            speaker_number=speaker_number+1;
           prompt={'请输入您的姓名'};
   name='输入姓名 ';
   numlines=1;
   defaultanswer={'无'};
   answer=inputdlg(prompt,name,numlines,defaultanswer);
   data(speaker_number).name=answer{1,1};                                   %保存用户信息
          data(speaker_number).means=nim1;
          data(speaker_number).cov=readcov(nis1);
          data(speaker_number).prob=nip1;
          data(speaker_number).pitch=s1;
          save('speech_database.dat','data','speaker_number','-append');
       else
           speaker_number=1;   
    prompt={'请输入用户姓名：'};
   name='用户 ';
   numlines=1;
   defaultanswer={'无'};
   answer=inputdlg(prompt,name,numlines,defaultanswer);                     %获得输入的名字
   data(speaker_number).name=answer{1,1};                                   %将输入的名字存入data数据表中
           data(speaker_number).means=nim1;
        data(speaker_number).cov=readcov(nis1);
        data(speaker_number).prob=nip1;
        data(speaker_number).pitch=s1;
           save('speech_database.dat','data','speaker_number');
       end
        
        message=strcat(answer{1,1},'的信息已成功存入数据库中 ')
        msgbox(message,'speechsignal DataBase','help')
    end
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %第二选项：声纹匹配
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
 fr=frm(strcat(pathname,namefile),16,8000,3);                               %读取音频进行降噪
 l=length(fr(1,:));                                                         %获取第一行的长度
 nosp=length(data);
 k=0;
 b=0;
 r=nosp;
 while(r~=1)                                                                %用k计算data中的数据量是2的多少次方
  r=floor(r/2);
  k=k+1;
end
p(2,nosp)=0;p(1,1)=0;
for i=1:nosp
    p(2,i)=i;
end
mc4=train(fr,filt,20);                                                      %提取梅尔频率倒谱系数，组成瀑布图
mc4=mc4(3:18,:);
mc=banshengsin(mc4);                                                        %为输入系数乘半升正弦函数
pitch2=pitch(pathname,namefile);                                            %求取基因周期
a=length(pitch2);
b=length(mc(1,:));
if a>b                                                                      %比较基音周期和说话声音的长度，基因周期长，那么去除基音周期多出部分的内容
    pitch2(b+1:a)=[];
else
    pitch2(a+1:b)=0;                                                        %如果说话声音长度长，为基音周期缺少的长度补0
end
mc=[mc;pitch2];                                                             %将基因周期补在最后一行
coff=length(mc(:,1));                                                       %获取第一列长度
o=length(mc(1,:));                                                          %获取第一行的长度
frameparts=struct('frame',{});
s=mod(l,k);
y=floor(l/k);
if s==0
   for i=1:k                                                                %data的长度是2的k次方
       frameparts(i).frame(coff,y)=0;                                       %以输入的音频的行数为行数，y为列数，给最后一个值复0
   end                                                                      %主要还是为了生成frameparts的矩阵把
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
c=length(data);                                                             %数据报表中成员的个数
for  i=1:k
   % tic
   p1=ident2(frameparts(i).frame,filt,data,p);                              %利用最大似然概率寻找输入数据的高斯模型
 %  toc
   p=upd_pr(p,p1);                                                          %迭代，用p1去更新p
   p=nmax1(p);                                                              %取p的最大值
end
p2=p(1)/o;
scores=zeros(nosp,1);
for i=1:nosp                                                                %nosp表中数据的长度
   pitch1=data(i).pitch';
  % tic
   scores(i,1)=myDTW(pitch2,pitch1(1:length(pitch2)));                      %返回表中数据与输入数据的最短路径
  % toc
end
scores;
[m,n]=sort(scores);                                                         %对scores升序排列

b=p(2,1);
if or((p2>-25),b==n)
nm=data(b).name;
       message=strcat('说话人是： ',nm);
       msgbox(message,'DataBase Info','help');
else
    message='无法匹配到相关人员';
    msgbox(message,'DataBase Info','help');
end
        else
            message='数据表为空，请添加用户';
            msgbox(message,'speech DataBase Error','warn');    
        end
        
    end 
   
%第三选项：删除全部数据，或只删除一个人的数据    
    if chos==3                                                  
        clc;
        close all;
        if (exist('speech_database.dat')==2)                                %当表中保存过任何人的声音就会返回2
             load('speech_database.dat','-mat');                            %将数据表用matlab标准格式打开
            button = questdlg('请问您想删除哪一位的数据?',...
                   '选择',...
                   '所有人','指定的一位','all');
             if strcmp(button,'所有人')                                     %判断按了哪个按钮
                delete('speech_database.dat');
                msgbox('您已成功删除所有数据','Database removed','help');
             else 
                  prompt={'请输入您想删除数据的人员姓名'};
   name='删除指定人员';
   numlines=1;
   defaultanswer={'无'};
   answer=inputdlg(prompt,name,numlines,defaultanswer);                     %生成输入框，answer为输入的名字
nspeaker=length(data)                   ;
 names=cell(1,nspeaker);
 for i=1:nspeaker                                                           %生成一个用于检索的向量
names{1,i}=data(i).name;
 end
[a,b]=ismember(answer{1,1},names);                                          %判断输入的名字是否在表中                            
 if a==0
       warndlg('不存在指定人选','警告！')
   else 
      data(b)=[];                                                           %找到了该名字人员所在位置，清空该处的音频
      speaker_number=length(data);
      save('speech_database.dat','data','speaker_number','-append');
   message=strcat('成功删除 ',answer{1,1},'的信息');
            msgbox(message,'成功删除','help');
 end
             end
        else
            warndlg('数据表中为空.',' 警告！ ')
        end
    end 
    %第四选项：音频信号分析
    if chos==4
        clc;
        close all;
        selezionato=0;
        while selezionato==0                                                %没有选择文件就一直选择
            [namefile,pathname]=uigetfile({'*.wav','speech signal (*.wav)'},'请选择音频文件');
            if namefile~=0
               [x,fs]=audioread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('请选择音频文件');
            end
        
        end
        %语谱图（瀑布图）
        figure('Name','瀑布图');
        Y1 = x(:,1);        %Y为双声道数据，取第2通道
        plot(Y1);           %画Y1波形图
        grid on;
        specgram(Y1,2048,44100,2048,1536);
         %Y1为波形数据
         %FFT帧长2048点(在44100Hz频率时约为46ms)
         %采样频率44.1KHz
         %加窗长度，一般与帧长相等
         %帧重叠长度，此处取为帧长的3/4
         T=1/fs;
         N=length(x);
         n=0:N-1;
         t=n/fs;
         figure('Name','频谱');
         subplot(2,2,1);plot(t,x);axis([0 N/fs -0.4 0.4]);
         title('时域信号');
         ylabel('振幅');xlabel('t');
         N0=256;
         n0=0:N0-1;
         f=n*fs/N;
         z=fft(x,N);mf1=abs(z');
         subplot(2,2,2);plot(f,mf1);xlabel('f');ylabel('幅值');
         subplot(2,2,3);plot(f,angle(z));xlabel('相位');ylabel('幅值');
         subplot(2,2,4);plot(freqz(x));


        
   %[im is ip]=init(mc1,16);
    % [nim1 nis1 nip1 times]=gmm(im,is,ip,mc1);
     
    end  
         if chos==5
        clc;
        close all;
        selezionato=0;
        while selezionato==0                                                %没有选择文件就一直选择
            [namefile,pathname]=uigetfile({'*.wav','speech signal (*.wav)'},'请选择音频文件');
            if namefile~=0
               [x,fs]=audioread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('请选择音频文件');
            end
        
        end
            x=x(:,1);
            audiowrite(namefile,x,fs);
            
         end
         %第六部分：语音去噪
              if chos==6
        clc;
        close all;
        selezionato=0;
        while selezionato==0                                                %没有选择文件就一直选择
            [namefile,pathname]=uigetfile({'*.wav','speech signal (*.wav)'},'请选择音频文件');
            if namefile~=0
               [sound,fs]=audioread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('请选择音频文件');
            end
        
        end
            sound=sound(:,1);                                               %读入单声道音频信号
            count=length(sound);
            t=(0:count-1)/fs;                                               %分段处理
            sound=sound';
            [c,l]=wavedec(sound,3,'db6');                                   %用小波函数db6对信号进行三层分解
            sigma=wnoisest(c,l,1);                                          %估计尺度为1的噪声标准偏差
            alpha=2;                                                        %调整参数，一般取2
            thr=wbmpen(c,l,sigma,alpha);                                    %获取消噪过程中的阈值
            keepapp=1;  
            yd=wdencmp('gbl',c,l,'db6',3,thr,'s',keepapp);                  %对信号进行消噪
            subplot(1,2,1)
            plot(t,sound);
            title('原始语音信号');
            axis([0 11 min(sound)-0.1 max(sound)+0.1])
            subplot(1,2,2)
            plot(t,yd);
            title('去噪语音信号');
            axis([0 11 min(yd)-0.1 max(yd)+0.1])
            audiowrite('final.wav',yd,fs);
              end
              %第七部分：声源分离
              if chos==7
        clc;
        close all;
        selezionato=0;
        while selezionato==0                                                %没有选择文件就一直选择,
            [namefile,pathname]=uigetfile({'*.wav','speech signal (*.wav)'},'请选择第一段音频文件');
            if namefile~=0
               [I1,fs]=audioread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('请选择音频文件');
            end
        
        end
        selezionato=0;  
        while selezionato==0                                                %没有选择文件就一直选择
            [namefile,pathname]=uigetfile({'*.wav','speech signal (*.wav)'},'请选择第二段音频文件');
            if namefile~=0
               [I2,fs]=audioread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('请选择音频文件');
            end
        
        end
        selezionato=0;  
        while selezionato==0                                                %没有选择文件就一直选择
            [namefile,pathname]=uigetfile({'*.wav','speech signal (*.wav)'},'请选择第四段音频文件');
            if namefile~=0
               [I3,fs]=audioread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('请选择音频文件');
            end
        
        end
        subplot(4,3,1),plot(I1),title('输入信号1');
        subplot(4,3,2),plot(I2),title('输入信号2');
        subplot(4,3,3),plot(I3),title('输入信号3');


        l1=length(I1);
        l2=length(I2);
        l3=length(I3);
        
        if l1>=l2&&l1>=l3                                                %维数统一化
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
                                                                         %将混合矩阵重新排列并输出
        subplot(4,3,4),plot(MixedS(1,:)),title('混合信号1');
        subplot(4,3,5),plot(MixedS(2,:)),title('混合信号2');
        subplot(4,3,6),plot(MixedS(3,:)),title('混合信号3');
       
        audiowrite('1mixwav1.wav',MixedS(1,:),fs);                       %保存wav数据
        audiowrite('1mixwav2.wav',MixedS(2,:),fs);
        audiowrite('1mixwav3.wav',MixedS(3,:),fs);
        MixedS_bak=MixedS;
%%%%%%%%%%%%%%%%%%%%%%%%%%标准化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%去均值，使矩阵为零均值
        MixedS_mean=zeros(3,1);
        for i=1:3                                                        %计算MixedS的均值
            MixedS_mean(i)=mean(MixedS(i,:));
        end                                                      

        for i=1:3                                                        %对矩阵中每一个值减去均值，使得矩阵变为零均值
             for j=1:size(MixedS,2)
                MixedS(i,j)=MixedS(i,j)-MixedS_mean(i);
             end
        end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%白化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%去除相关性，可以降低问题的复杂性，使得算法收敛性变好，简化后续独立分量的提取问题
            
        MixedS_cov=cov(MixedS');                                         %cov为求协方差的函数
        [E,D]=eig(MixedS_cov);                                           %对信号矩阵的协方差函数进行特征值分解
        %Q=inv(sqrt(D))*(E)';                                            %Q为白化矩阵
        Q=sqrt(D)\(E)';
        MixedS_white=Q*MixedS;                                           %MixedS_white为白化后的信号矩阵
        IsI=cov(MixedS_white');                                          %IsI应为单位阵
	
%%%%%%%%%%%%%%%%%%%%%%%%FASTICA算法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%本质是寻找声源的非高斯最大的不动点迭代方案
        X=MixedS_white;                                                  %以下算法将对X进行操作
        [VariableNum,SampleNum]=size(X);                                 
        numofIC=VariableNum;                                             %在此应用中，独立元个数等于变量个数(即声源个数)
        B=zeros(numofIC,VariableNum);                                    %初始化列向量w的寄存矩阵,B=[b1 b2 ... bd]
        for	r=1:numofIC                                                  %根据声源的个数算出相应个数的独立声源
            i=1;maxIterationsNum=100;                                    %设置最大迭代次数（即对于每个独立分量而言迭代均不超过此次数）
            IterationsNum=0;	
            b=rand(numofIC,1)-.5;                                        %随机设置b初值,
            b=b/norm(b);                                                 %对b标准化norm(b):向量元素平方和开根号,对b进行归一化
            while i<=maxIterationsNum+1                                  %进行多次迭代，逼近接收器的逆矩阵
                if i == maxIterationsNum                                 %若在最大次数内不收敛，循环结束处理
                        fprintf('\n第%d分量在%d次迭代内并不收敛。', r,maxIterationsNum); 
                        break	;
                end	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%核心公式%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                bOld=b;                                 
                a2=1;	
                u=1;	
                t=X'*b;	
                g=t.*exp(-a2*t.^2/2);	
                dg=(1-a2*t.^2).*exp(-a2*t.^2/2);	
                b=((1-u)*t'*g*b+u*X*g)/SampleNum-mean(dg)*b;	
                b=b-B*B'*b;                                              %对b正交化
                b=b/norm(b);                                             %对b归一化
                if abs(abs(b'*bOld)-1)<1e-9                              %如果足够近似，则
                        B(:,r)=b;                                        %保存所得向量b
                        break;
                end 
            i=i+1; 
             end	
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%ICA计算的数据复原并构图%%%%%%%%%%%%%%%%%%%%%%%%%
    ICAedS=B'*Q*MixedS_bak;                                         %计算ICA后的矩阵，得到的信号求逆对源信号进行分解
                                                                    %将混合矩阵重新排列并输出
    subplot(4,3,7),plot(ICAedS(1,:)),title('ICA解混信号1');
    subplot(4,3,8),plot(ICAedS(2,:)),title('ICA解混信号2');
    subplot(4,3,9),plot(ICAedS(3,:)),title('ICA解混信号3');

    audiowrite('1dstwav1.wav',ICAedS(1,:),fs);                      %保存wav数据
    audiowrite('1dstwav2.wav',ICAedS(2,:),fs);
    audiowrite('1dstwav3.wav',ICAedS(3,:),fs);

              end
end