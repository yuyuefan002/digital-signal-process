
%����һmelƵ���˲����飬ÿ���˲�������һ�������˲���
%spacing���˲������ͨ�����
%bandwidth���˲������ͨ��
%p�Ǵ������˲����ĸ���
%����һ��n*p����ÿһ�д���һ���˲�������Ӧ
function filt=melfilter(spacing,bandwidth,p)
mcentre=bandwidth/2;
centre=melinv(mcentre);
lastf=melinv(p*bandwidth);
for i=1:p
	mstart=mcentre-bandwidth/2;
	mend=mcentre+bandwidth/2;
	start=melinv(mstart);
	term=melinv(mend);
	trfilt=triang(round(term)-round(start)+1);
	c=0;
	for j=round(start):round(term)
		filt(j+1,i)=trfilt(c+1);
		c=c+1;
	end
	mcentre=mcentre+spacing;
	centre=melinv(mcentre);
end



	
