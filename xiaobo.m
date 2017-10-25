

clear all
close all

%������ǿ����
%����Ϊ���Ը�˹������

%��������������sound�С�
[sound,fs]=audioread('sound.wav');
sound=sound(:,1);
count=length(sound);
t=(0:count-1)/fs;

%��������ļ��봦��
noise=0.05*randn(1,count);
y=sound'+noise;
wavwrite(y,fs,'compound.wav');

%��С������db6���źŽ�������ֽ�
[c,l]=wavedec(y,3,'db6');

%���Ƴ߶�Ϊ1��������׼ƫ��
sigma=wnoisest(c,l,1);
alpha=2;

%��ȡ��������е���ֵ
thr=wbmpen(c,l,sigma,alpha);
keepapp=1;

%���źŽ�������
yd=wdencmp('gbl',c,l,'db6',3,thr,'s',keepapp);
subplot(1,2,1)
plot(t,sound);
title('ԭʼ�����ź�');
axis([0 11 min(sound)-0.1 max(sound)+0.1])

subplot(1,2,2)
plot(t,yd);
title('ȥ�������ź�');
axis([0 11 min(yd)-0.1 max(yd)+0.1])

wavwrite(yd,fs,'final.wav');

