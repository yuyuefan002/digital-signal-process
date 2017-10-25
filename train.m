function mels=train(snd,filt,nmfcc)
n=length(snd(1,:));                                                         %获取一维音频信号长度
for i=1:n
%	fft_snd=fft(snd(:,i),10423);
%	fft_snd=fft_snd(1:4000);
%	fft_snd(5211)=0;
	fft_snd=fft(snd(:,i),10423);                                            %对每一帧做fft
	fft_snd(1)=[];                                                          %不需要直流分量
	fft_snd=fft_snd(1:5211);                                                %只取前一半频率的值,因为频域图在其中点对称
	pwr_snd=abs(fft_snd).^2;                                                %获取功率谱
	tmels=mfcc(filt,pwr_snd,nmfcc);                                         %获取功率谱的包络，此为梅尔频率倒谱系数
	mels(:,i)=tmels';                                                       %将其转置赋保存，这样每一列就是一个时段的瀑布图
end
	