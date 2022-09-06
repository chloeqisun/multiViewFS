function ie = kInformationEntopy(sig,SegmentNum)
% 计算信号的信息熵
% 参考《矿用带式输送机托辊远程故障诊断系统》
% 输入：
% sig：输入信号
% SegmentNum：拟分组数，如果不输入，则自动使用斯特格斯（Sturges）经验公式计算。
% 输出：
% ie：信息熵求解结果
[len,num] = size(sig);
if num == 2
    SigLen = length(sig);  %输入信号长度
    if nargin == 1
        SegmentNum = round(1.87*(SigLen-1)^(2/5));  %斯特格斯（Sturges）的经验公式，求最佳分组数
    end
    CutLen = SigLen/SegmentNum; %每组信号长度
    Ent = [];
    for i = 1:SegmentNum
        Ent = [Ent;sum(sig(round(CutLen*(i-1)+1):round(CutLen*i),:),1)];  %求每组信号的总能量
    end
    pk = Ent./sum(Ent,2); %第k段能量占总能量大小
    ie = -sum(pk.*log(pk)); %信息熵公式
end 
if num == 1
    SigLen = length(sig);  %输入信号长度
    if nargin == 1
        SegmentNum = round(1.87*(SigLen-1)^(2/5));  %斯特格斯（Sturges）的经验公式，求最佳分组数
    end
    CutLen = SigLen/SegmentNum; %每组信号长度
    for i = 1:SegmentNum
        Ent(i) = sum(sig(round(CutLen*(i-1)+1):round(CutLen*i)));  %求每组信号的总能量
    end
    pk = Ent/sum(Ent); %第k段能量占总能量大小
    ie = -sum(pk.*log(pk)); %信息熵公式
end
end

