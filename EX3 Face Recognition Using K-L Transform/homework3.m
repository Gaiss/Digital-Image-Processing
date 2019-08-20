clear all
%% 读取数据
A = [];  % 训练集
A_label = [];  % 训练集标签
B = [];  % 测试集
B_label = [];  % 测试集标签
for i = 1:15
    for j = 1:8 % 每人取8张做训练
        num = 10000 + i*100 + j;
        num = num2str(num);
        num(1) = '0';
        f = imread(strcat('yale/',num,'.bmp'));
        f = imresize(f,[64 64]);  % 压缩图片，提高运行速度
        A = [A,f(:)];
        A_label = [A_label,i];
    end
    for j = 9:11 % 每人取3张做测试
        num = 10000 + i*100 + j;
        num = num2str(num);
        num(1) = '0';
        f = imread(strcat('yale/',num,'.bmp'));
        f = imresize(f,[64 64]);
        B = [B,f(:)];
        B_label = [B_label,i];
    end
end
A = double(A); % 训练集，每列为一个样本
B = double(B); % 测试集，每列为一个样本
%% K-L变换 + PCA降维
C = cov(A'); % 总体散布矩阵
[V D] = eig(C);
D = diag(D); % 特征值
[Ds index] = sort(D,'descend'); % 降序排列特征值
for k = 1:length(Ds) % PCA 保留前p个特征值占总能量的92%
    if sum(Ds(1:k))/sum(Ds) > 0.92
        p = k;
        break
    end
end
Vs = V(:,index);
W = Vs(:,1:p); % 投影矩阵
Y_A = W'*A;  % 训练集特征提取
Y_B = W'*B;  % 测试集特征提取
%% 相似度计算
[~,I1] = pdist2(Y_A',Y_B','seuclidean','Smallest',1);  % 法一，标准化欧氏距离,取最小
[~,I2] = pdist2(Y_A',Y_B','cosine','Smallest',1);  % 法二，余弦距离，取最小
B_testlabel1 = A_label(I1); % 法一预测测试集标签
B_testlabel2 = A_label(I2); % 法二
Accuracy1 = sum(B_testlabel1 == B_label)/length(B_label) % 法一准确率
Accuracy2 = sum(B_testlabel2 == B_label)/length(B_label) % 法二准确率
%% 输出具体识别结果
NUM = [];
for i = 1:15
    for j = 9:11
        num = 10000 + i*100 + j;
        num = int2str(num);
        num(1) = '0';
        NUM = [NUM;num];
    end
end
iftrue1 = (B_testlabel1 == B_label);
iftrue2 = (B_testlabel2 == B_label);
table = table(NUM,B_label',B_testlabel1',iftrue1',B_testlabel2',iftrue2',...
    'VariableNames',{'Name','real','predict1','iftrue1','predict2','iftrue2'})