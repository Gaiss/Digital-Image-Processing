clear all
%% ��ȡ����
A = [];  % ѵ����
A_label = [];  % ѵ������ǩ
B = [];  % ���Լ�
B_label = [];  % ���Լ���ǩ
for i = 1:15
    for j = 1:8 % ÿ��ȡ8����ѵ��
        num = 10000 + i*100 + j;
        num = num2str(num);
        num(1) = '0';
        f = imread(strcat('yale/',num,'.bmp'));
        f = imresize(f,[64 64]);  % ѹ��ͼƬ����������ٶ�
        A = [A,f(:)];
        A_label = [A_label,i];
    end
    for j = 9:11 % ÿ��ȡ3��������
        num = 10000 + i*100 + j;
        num = num2str(num);
        num(1) = '0';
        f = imread(strcat('yale/',num,'.bmp'));
        f = imresize(f,[64 64]);
        B = [B,f(:)];
        B_label = [B_label,i];
    end
end
A = double(A); % ѵ������ÿ��Ϊһ������
B = double(B); % ���Լ���ÿ��Ϊһ������
%% K-L�任 + PCA��ά
C = cov(A'); % ����ɢ������
[V D] = eig(C);
D = diag(D); % ����ֵ
[Ds index] = sort(D,'descend'); % ������������ֵ
for k = 1:length(Ds) % PCA ����ǰp������ֵռ��������92%
    if sum(Ds(1:k))/sum(Ds) > 0.92
        p = k;
        break
    end
end
Vs = V(:,index);
W = Vs(:,1:p); % ͶӰ����
Y_A = W'*A;  % ѵ����������ȡ
Y_B = W'*B;  % ���Լ�������ȡ
%% ���ƶȼ���
[~,I1] = pdist2(Y_A',Y_B','seuclidean','Smallest',1);  % ��һ����׼��ŷ�Ͼ���,ȡ��С
[~,I2] = pdist2(Y_A',Y_B','cosine','Smallest',1);  % ���������Ҿ��룬ȡ��С
B_testlabel1 = A_label(I1); % ��һԤ����Լ���ǩ
B_testlabel2 = A_label(I2); % ����
Accuracy1 = sum(B_testlabel1 == B_label)/length(B_label) % ��һ׼ȷ��
Accuracy2 = sum(B_testlabel2 == B_label)/length(B_label) % ����׼ȷ��
%% �������ʶ����
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