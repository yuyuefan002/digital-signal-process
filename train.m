function mels=train(snd,filt,nmfcc)
n=length(snd(1,:));                                                         %��ȡһά��Ƶ�źų���
for i=1:n
%	fft_snd=fft(snd(:,i),10423);
%	fft_snd=fft_snd(1:4000);
%	fft_snd(5211)=0;
	fft_snd=fft(snd(:,i),10423);                                            %��ÿһ֡��fft
	fft_snd(1)=[];                                                          %����Ҫֱ������
	fft_snd=fft_snd(1:5211);                                                %ֻȡǰһ��Ƶ�ʵ�ֵ,��ΪƵ��ͼ�����е�Գ�
	pwr_snd=abs(fft_snd).^2;                                                %��ȡ������
	tmels=mfcc(filt,pwr_snd,nmfcc);                                         %��ȡ�����׵İ��磬��Ϊ÷��Ƶ�ʵ���ϵ��
	mels(:,i)=tmels';                                                       %����ת�ø����棬����ÿһ�о���һ��ʱ�ε��ٲ�ͼ
end
	