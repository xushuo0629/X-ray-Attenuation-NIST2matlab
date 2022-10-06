注意：./data 文件夹 必须与 CalculateMassAC.mexw64同一个目录

1  调用方法:
	R = CalculateMassAC( formula,  E )
	formula 是指定元素组成的char字符矩阵表示的化学式 或者 2行N列数值矩阵

1.1 输入参数 formula 
formula是必选参数
如果formula是字符矩阵，则需要满足：
原子总是以一个大写字母开始，接着跟随 0 个或任意个小写字母，表示原子的名字，比如: Be, Pb, C.
如果数量大于 1，原子后会跟着数字表示原子的数量。如果数量等于 1 则省略或者不省略都可以。
例如，"H2O" 和 "H2O2" 是可行的，"H1O2" 这个表达也是可行的。
两个化学式连在一起可以构成新的化学式。例如 "H2O2He3Mg4" 也是化学式。
由括号括起的化学式并佐以数字（可选择性添加）也是合法化学式。例如 "(H2O2)" 和 "(H2O2)3" 是化学式
数字可以是小数形式，比如"H1.0O2.0"、"C0.056H0.126K.2"是合法的，但不能是科学计数法形式，比如"H1e3O1e-2"是不合法的
用这个规则可以表示混合物例如 碳和硅的2：5混合物可以表示为"C2Si5"或者"C20Si50"都可以，因为采用独立原子模型时
质量衰减系数几乎只与元素配比有关，与化学键怎么连的无关，当考虑比较仔细时，比如考虑原子之间的化学键对质量衰减系数的影响时
本程序不再适用。
在MATLAB里字符矩阵必须写成  'H2O' 即单引号形式，不能写成双引号"H2O"

如果formula 是数值矩阵，则，必须是2*N的，第一行对应原子序数，第二行对应元素的原子数比例 
例如: formula = [1, 8 ; 2, 1] 代表“H2O”,  [1, 8 ; 0.2, 0.1]也代表“H2O”,因为只关心相对比例，不关心决定数值

1.2 能量E
E 是可选参数 如果没有E，则采用默认的能量点

E必须是形如1*N的数值矩阵，即一行N列，单位必须是keV

1.3 输出参数 R
一个结构体 ，数组每个域名代表一种衰减系数，单位都是cm^2/g
要得到线性衰减系数，只需要自己查询相关材料密度即可。
域名含义：
PhotoelectricEdge: K-edge L-edge等对应的原子序数
PhotonEnergy： 光子能量/keV
CoherentScattering: 相干散射
IncoherentScattering: 非相干散射，即康普顿散射
PhotoelectricAbsorption: 光电效应
NuclearFieldPairProduction：原子核电场里的电子对效应
ElectronFieldPairProduction：电子电场所产生的电子对效应
TotalWithCoherent：含相干散射的总散射
TotalWithoutCoherent: 不含相干散射的总散射

调用示例(化学式是示意，未必真实存在相应物质)

R1 = CalculateMassAC('Pb2O3(H2O).7(CuO)6.5' , 10:0.01:200)

'Pb2O3(H2O).7(CuS)6.5' 表示该物质中含有2个Pb，3个O，0.7个H2O，6.5个CuS
经过合并，可以等价于输入 'Pb2O3.7H1.4Cu6.5S6.5'
10:0.01:200 表示10keV到200keV，间隔为0.01keV
能量最低为1keV，最高到GeV, 在X射线领域已经足够了

2 Setup
如果修改了源文件，需要重新编译
需要将CalculateMassAC.cpp所在目录设置为MATLAB当前目录，并在MATLAB命令行输入

mex -setup c++
mex CalculateMassAC.cpp xcom.cpp