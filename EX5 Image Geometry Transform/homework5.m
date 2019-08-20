%% 法一 warp函数
clear all
A = imread('Colorful Rose.jpg');
xx = 0:0.01:10;
yy = -10:0.01:10;
[x1 y1] = meshgrid(xx,yy);
z1 = x1.^2 + y1.^2;
z1(z1<10) = NaN;
z1(z1>100) = NaN;
% 
% figure()
% subplot(1,2,1),imshow(A)
% title('Original')
% subplot(1,2,2),warp(x,y,z,A);
% view(3)
% axis square
% grid on
% title('warp函数')
% 
% figure()
% subplot(2,2,1),warp(x,y,z,A)
% title('主视图')
% view(0,0)
% axis square
% subplot(2,2,2),warp(x,y,z,A)
% title('左视图')
% view(90,0)
% axis square
% subplot(2,2,3),warp(x,y,z,A)
% title('俯视图')
% view(2)
% axis square
%% 法二 CData属性 前向映射
% clear all
% A = imread('Colorful Rose.jpg');
% [m n ~] = size(A);
% xx = linspace(0,10,n);
% yy = linspace(10,-10,m);
% [x2 y2] = meshgrid(xx,yy);
% z2 = x2.^2 + y2.^2;
% B = z2<=100 & z2>=10;
% z2(B == 0) = NaN;
% 
% figure()
% subplot(1,2,1),imshow(A)
% title('Original')
% subplot(1,2,2)
% h=surf(x1,y1,z1);
% shading interp
% set(h,'CData',A); % 纹理贴图
% view(3)
% axis square
% title('CData属性 前向映射')
% 
% figure()
% subplot(2,2,1)
% h=surf(x1,y1,z1);
% shading interp
% set(h,'CData',A);
% title('主视图')
% view(0,0)
% axis square
% subplot(2,2,2)
% h=surf(x1,y1,z1);
% shading interp
% set(h,'CData',A);
% title('左视图')
% view(90,0)
% axis square
% subplot(2,2,3)
% h=surf(x1,y1,z1);
% shading interp
% set(h,'CData',A);
% title('俯视图')
% view(2)
% axis square
%% 法三 CData属性 后向映射
% clear all
A = imread('Colorful Rose.jpg');
[m n ~] = size(A);
theta = linspace(-pi/2,pi/2,110) ;
r = linspace(sqrt(10),10,73);
[theta r] = meshgrid(theta,r);
x = r.*cos(theta) ;
y = r.*sin(theta) ;
z = x.^2 + y.^2;

% figure()
% subplot(1,2,1),imshow(A)
% title('Original')
% subplot(1,2,2)
% mesh(x,y,z,'EdgeColor','k')
% axis off
% shading interp
% set(h,'CData',A); % 纹理贴图
% view(3)
% axis square
% title('CData属性 后向映射')
% 
% figure()
% subplot(2,2,1)
% h=surf(x,y,z);
% shading interp
% set(h,'CData',A);
% title('主视图')
% view(0,0)
% axis square
% subplot(2,2,2)
% h=surf(x,y,z);
% shading interp
% set(h,'CData',A);
% title('左视图')
% view(90,0)
% axis square
% subplot(2,2,3)
% h=surf(x,y,z);
% shading interp
% set(h,'CData',A);
% title('俯视图')
% view(2)
% axis square
%%
% figure()
% subplot(1,3,1),warp(x1,y1,z1,A);
% view(2)
% axis square
% grid on
% title('warp函数')
% subplot(1,3,2)
% h=surf(x2,y2,z2);
% shading interp
% set(h,'CData',A); % 纹理贴图
% view(2)
% axis square
% title('CData属性 前向映射')
% subplot(1,3,3)
% h=surf(x,y,z);
% shading interp
% set(h,'CData',A); % 纹理贴图
% view(2)
% axis square
% title('CData属性 后向映射')

% w=5; %控制分块窗口的大小，16×16
% M=m/w;N=n/w;  %行列间隔块个数
% % 区域块分割
% t1=(0:M-1)*w+1;t2=(1:M)*w;
% t3=(0:N-1)*w+1;t4=(1:N)*w;
% figure() 
% imshow(A); 
% hold on
% for i=1:M
%     for j=1:N
%         x=t1(i):t2(i);
%         y=t3(j):t4(j);
%         rectangle ('Position', [t3(j) t1(i) length(x) length(y)], ...,
%             'EdgeColor','r','LineWidth',0.05); %绘制矩形分割格子，红色
%     end
% end