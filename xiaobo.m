

clear all
close all

%语音增强程序
%噪声为加性高斯白噪声

%读入声音到数组sound中。
[sound,fs]=audioread('sound.wav');
sound=sound(:,1);
count=length(sound);
t=(0:count-1)/fs;

%完成声音的加噪处理
noise=0.05*randn(1,count);
y=sound'+noise;
wavwrite(y,fs,'compound.wav');

%用小波函数db6对信号进行三层分解
[c,l]=wavedec(y,3,'db6');

%估计尺度为1的噪声标准偏差
sigma=wnoisest(c,l,1);
alpha=2;

%获取消噪过程中的阈值
thr=wbmpen(c,l,sigma,alpha);
keepapp=1;

%对信号进行消噪
yd=wdencmp('gbl',c,l,'db6',3,thr,'s',keepapp);
subplot(1,2,1)
plot(t,sound);
title('原始语音信号');
axis([0 11 min(sound)-0.1 max(sound)+0.1])

subplot(1,2,2)
plot(t,yd);
title('去噪语音信号');
axis([0 11 min(yd)-0.1 max(yd)+0.1])

wavwrite(yd,fs,'final.wav');

