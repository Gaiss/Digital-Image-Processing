% clear all
% A=imread('图片1.png');
% A=A(:,:,1);
% [m n]=size(A);
% B=zeros(m,n,3);
% for x=1:m
%     for y=1:n
%         if A(x,y)<=63
%             B(x,y,1)=0;
%             B(x,y,2)=254-4*A(x,y);
%             B(x,y,3)=255;
%         else if A(x,y)<=127
%                B(x,y,1)=0;
%                B(x,y,2)=4*A(x,y)-254;
%                B(x,y,3)=510-4*A(x,y); 
%         else if A(x,y)<=191
%                B(x,y,1)=4*A(x,y)-510;
%                B(x,y,2)=255;
%                B(x,y,3)=0; 
%             else 
%                 B(x,y,1)=255;
%                 B(x,y,2)=1022-4*A(x,y);
%                 B(x,y,3)=0;
%             end
%             end
%         end
%     end
% end
% figure()
% subplot(1,2,1),imshow(A)
% title('灰度图像')
% subplot(1,2,2),imshow(uint8(B));
% title('伪彩色处理')

%% 图像划线分割
clc; clear all; close all; 
A=imread('图片1.png');
A=A(:,:,1);
A=imresize(A,[128,160]);
[m n]=size(A);  
w=16; %控制分块窗口的大小，16×16
M=m/w;N=n/w;  %行列间隔块个数
% 区域块分割
t1=(0:M-1)*w+1;t2=(1:M)*w;
t3=(0:N-1)*w+1;t4=(1:N)*w;
figure() 
imshow(A); 
hold on
for i=1:M
    for j=1:N
        x=t1(i):t2(i);
        y=t3(j):t4(j);
        rectangle ('Position', [t3(j) t1(i) length(x) length(y)], ...,
            'EdgeColor','r','LineWidth',0.05); %绘制矩形分割格子，红色
    end
end
%% 图像分块
% A=imread('图片1.png');
% A=A(:,:,1);
% A=imresize(A,[128,160]);
% [m n]=size(A); 
% w=16; %窗口大小
% M=m/w;N=n/w;
% t1=(0:M-1)*w+1;t2=(1:M)*w;
% t3=(0:N-1)*w+1;t4=(1:N)*w;
% figure()
% k=0;
% for i=1:M
%     for j=1:N      
%         temp=A(t1(i):t2(i),t3(j):t4(j),:);
%         k=k+1;
%         subplot(M,N,k);
%         imshow(temp);
%     end
% end
%%
% A='lenna.jpg';
% B='Lenna_(test_image).png';
% R=gray2rgb(A,B);
% function R=gray2rgb(img1,img2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %　This Program converts a gray image ro RGB image based on the colors of the destination image. The　better the destination image match with the source gray image, the better the coloring will be. The program takes some time  as the searching time is high. You can decrease the searching time by taking　only samples from the used color image but quality may decrease. U can use jittered sampling for improving running speed. % 
% %  You can use also use the attahed test images, Use the following combinations for better result　nature1.jpg(as img1) and nature2.jpg(as img2) or test1.jpg(as img1) and test2.jpg (as img2)　%
% %　Usage: gray2rgb('nature1.jpg','nature2.jpg');  %
% 
% %  Authors : Jeny　Rajan , Chandrashekar　P.S %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % img1 - Source Image  (gray image)   
% % img2 - Selected color image for coloring the gray image. 
% tic
%  clc;
%  warning off;
%  imt=imread(img1);
%  ims=imread(img2);
%  [sx sy sz]=size(imt);
%  [tx ty tz]=size(ims);
%  if sz~=1
%      imt=rgb2gray(imt);
%  end
%  if tz~=3
%      disp ('img2 must be a color image (not indexed)');
%  else
%      imt(:,:,2)=imt(:,:,1);
%      imt(:,:,3)=imt(:,:,1);
% 
%  % Converting to ycbcr color space
%      nspace1=rgb2ycbcr(ims);
%      nspace2= rgb2ycbcr(imt);
% 
%     ms=double(nspace1(:,:,1));
%      mt=double(nspace2(:,:,1));
%      m1=max(max(ms));
%      m2=min(min(ms));
%      m3=max(max(mt));
%      m4=min(min(mt));
%      d1=m1-m2;
%      d2=m3-m4;
%  % Normalization
%      dx1=ms;
%      dx2=mt;
%      dx1=(dx1*255)/(255-d1);
%      dx2=(dx2*255)/(255-d2);
%      [mx,my,mz]=size(dx2);
%  %Luminance Comparison
%      disp('Please wait..................');
%      for i=1:mx
%          for j=1:my
%               iy=dx2(i,j);
%               tmp=abs(dx1-iy);
%               ck=min(min(tmp));
%               [r,c] = find(tmp==ck);
%               ck=isempty(r);
%               if (ck~=1)            
%                   nimage(i,j,2)=nspace1(r(1),c(1),2);
%                   nimage(i,j,3)=nspace1(r(1),c(1),3);
%                   nimage(i,j,1)=nspace2(i,j,1);           
%              end
%           end
%       end
%      rslt=ycbcr2rgb(nimage);
%      figure,imshow(uint8(imt));
%      figure,imshow(uint8(rslt));
%      R=uint8(rslt);
%      toc
%  end  
% end

% x=0:255
% for i=1:64
%     y1(i)=0;
%     y2(i)=254-4*x(i);
%     y3(i)=255;
% end
% for i=65:128
%         y1(i)=0;
%         y2(i)=4*x(i)-254;
%         y3(i)=510-4*x(i); 
% end
% for i=129:192          
%     y1(i)=4*x(i)-510;
%             y2(i)=255;
%             y3(i)=0; 
% end
% for i=193:256
%             y1(i)=255;
%             y2(i)=1022-4*x(i);
%             y3(i)=0;
%         end
% 
% figure()
% hold on
% plot(y1,'r-','linewidth',2)
% plot(y2,'g-','linewidth',2)
% plot(y3,'b-','linewidth',2)
% legend('R(x,y)','G(x,y)','B(x,y)')
% axis([0 270 0 320])
% xlabel('f(x,y)')
% ylabel('RGB')