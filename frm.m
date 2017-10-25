%预处理，去掉没有声音的部分，降噪，分帧
function fr=frm(file,duration,fs,x)
nsig=audioread(file);
sig=trim(nsig);
if length(sig)>400000
	sig=sig(1:400000);
end
sig=sig(1:round(length(sig)/x));                                            %音频信号分帧
[thr,sorh,keepapp] = ddencmp('den','wv',sig);                               %为使用小波去噪做准备，返回阈值
sig = wdencmp('gbl',sig,'db3',2,thr,sorh,keepapp);                          %小波去噪
l=length(sig);
width=fs*duration/1000;
nframes=floor((3*l/(2*width))-2);                                           %算一下有几帧
h=hamming(width-1);                                                         %设置汉明窗，使用汉明窗是为了截取一定范围的信号，同时又不造成频率泄露
check=1;
numf=1;
crossperfrm=14;
fr(width-1,1)=0;
for i=1:nframes
	t=floor((i-1)*2*width/3)+1;
	count=0;
	if (t+width-2)>l 
		for j=t:l
			if check*sig(j)<0
				count=count+1;
				check=sig(j)/abs(sig(j));
			end
		end
		if count<crossperfrm
			fr(:,numf)=sig(t:l);
			numf=numf+1;
		end
	else
		for j=t:t+width-2
			if check*sig(j)<0
				count=count+1;
				check=sig(j)/abs(sig(j));
			end
		end
		if count<crossperfrm
			fr(:,numf)=sig(t:t+width-2);
			numf=numf+1;
		end
	end
	if count<crossperfrm
		fr(:,numf-1)=fr(:,numf-1).*h;
	end
end