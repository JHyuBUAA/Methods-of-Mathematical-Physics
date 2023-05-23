IB = imread('Brain_CT.bmp');
IB = im2double(IB);
I0 = rgb2gray(IB);%读入灰度图像
[m,n]=size(I0);
I1=Butterworth_low(I0);
figure
II0=fftshift(fft2(I0));
II0=log(abs(II0)+1);
II1=fftshift(fft2(I1));%二维傅里叶变换
II1=log(abs(II1)+1);%获取幅度谱
subplot(221)
imshow(I0);
title("原图")
subplot(222)
imshow(I1);
title("低通滤波")
subplot(223)
imshow(II0);
title("原图幅度谱")
subplot(224)
imshow(II1)
title("滤波后幅度谱")
grow_res=grow(I1);
for xx=1:m
    for yy=1:n
        if grow_res(xx,yy)==1
            IB(xx,yy,:)=[1,0,0];
        end
    end
end
figure
subplot(121)
imshow(grow_res)
title("病灶区域")
subplot(122)
imshow(IB)
title("原图中的病灶区域")

%%区域生长法分割
function grow1 = grow(I1)
figure 
imshow(I1)
[M,N]=size(I1);
[y,x]=getpts; %单击取点后，按enter结束
x1=round(x);%获得整数x
y1=round(y);%获得整数y
seed=I1(x1,y1); %获取中心像素灰度值
 
grow1=zeros(M,N);%要绘制二值图，就可以先设置为0数组
grow1(x1,y1)=1;%种子点设为1，满足条件的点在后续判断中也会设为一
 
growcount=1; %待处理点个数
threshold=0.15;%与isinteger相互结合起来理解，如果没有isinteger，该值应该变得很大
 
while growcount>0
    growcount=0;
     for i=1:M %让种子去融合图片中满足的所有点
        for j=1:N
            if grow1(i,j)==1 %点在"栈"内，即一开始先检测到种子点，然后根据判断条件去判断种子点邻域的点是否满足条件，若满足为1，然后进行下一轮循环进行判断，直到循环时没有点可以满足添加进行融合了
                if (i-1)>1&(i+1)<M&(j-1)>1&(j+1)<N %确保可以将8邻域完全放置
                    for u=-1:1 %8邻域生长
                        for v=-1:1
                            if (grow1(i+u,j+v)==0&abs(I1(i+u,j+v)-seed)<=threshold)%只有在点为0的时候才进行，否则gowtcount会计算错误
                                grow1(i+u,j+v)=1;%种子合并其他点，只要满足条件的点都会和种子点一样置1
                                growcount=growcount+1;  %记录此次新生长的点个数
                            end
                        end
                    end
                end
            end
        end
    end
end
 
end


%%巴特沃斯低通滤波器
function I = Butterworth_low(I_in)
    f=double(I_in);
    g=fft2(f);     %采用傅立叶变换
    g=fftshift(g); %数据矩阵平移，将变换后的图像频谱中心从矩阵的原点移到矩阵的中心
    [N1,N2]=size(g);
    result=zeros(N1,N2);%将结果存在result里
    n=2;
    d0=60;
    n1=fix(N1/2);
    n2=fix(N2/2);
    H=zeros(N1,N2);
    for i=1:N1
        for j=1:N2
            d=sqrt((i-n1)^2+(j-n2)^2);%频率点与频域中心的距离
            h=1/(1+0.414*(d/d0)^(2*n));%计算Butterworth低通变换函数
            H(i,j)=h;%记录滤波器
            result(i,j)=h*g(i,j);%频域相乘获得结果
        end
    end
    result=ifftshift(result);
    X2=ifft2(result);%反变换获得结果
    I=real(X2);
    figure
    imshow(H)%将滤波器显示为图像
    title("滤波器")
end


