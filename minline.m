function datamin=minline(data,linenum)
len=size(data,1);
datamin=zeros(len,1);
for i=1:len
    y=data(i,:);
    datamin(i)=min(y(end-linenum(i)+1:end));
end
