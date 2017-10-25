%Ԥ����ȥ��û�������Ĳ��֣����룬��֡
function fr=frm(file,duration,fs,x)
nsig=audioread(file);
sig=trim(nsig);
if length(sig)>400000
	sig=sig(1:400000);
end
sig=sig(1:round(length(sig)/x));                                            %��Ƶ�źŷ�֡
[thr,sorh,keepapp] = ddencmp('den','wv',sig);                               %Ϊʹ��С��ȥ����׼����������ֵ
sig = wdencmp('gbl',sig,'db3',2,thr,sorh,keepapp);                          %С��ȥ��
l=length(sig);
width=fs*duration/1000;
nframes=floor((3*l/(2*width))-2);                                           %��һ���м�֡
h=hamming(width-1);                                                         %���ú�������ʹ�ú�������Ϊ�˽�ȡһ����Χ���źţ�ͬʱ�ֲ����Ƶ��й¶
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