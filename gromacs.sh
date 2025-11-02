2025/10/29
###DNA高度标准化的生物大分子构件
#GROMACS 分子动力学模拟
#去水分子
grep -v HOH 1s40.pdb > 1s40_clean.pdb
#top结构
gmx pdb2gmx -f 1s40_clean.pdb -o 1s40_clean.gro -p 1s40_clean.top -i 1s40_clean.itp  -water tip3p -ignh##忽略氢键，选3，94力场
##添加盒子
gmx editconf -f 1s40_clean.gro -o 1s40_box.gro -c -d 1.0 -bt dodecahedron
gmx editconf -f 1s40.gro -o 1s40box.gro -c -d 1.0 -bt dodecahedron
###填充溶剂
gmx solvate -cp 1s40_box.gro -cs spc216.gro -o 1s40_solv.gro -p 1s40_clean.top
gmx solvate -cp 1s40box.gro -cs spc216.gro -o 1s40solv.gro -p 1s40.top
###离子嵌入
gmx grompp -f ions.mdp -c 1s40_solv.gro -p 1s40_clean.top -o 1s40_ions.tpr -maxwarn 1
gmx grompp -f ions.mdp -c 1s40solv.gro -p 1s40.top -o 1s40ions.tpr -maxwarn 1
##替换溶剂为离子
gmx genion -s 1s40ions.tpr -o 1s40_solv_ions.gro -p 1s40.top -pname NA -nname CL -neutral##14 SOL
###能量最小化：创建em.mdp
; 运行控制
integrator              = steep       ; 最速下降法（适用于初始EM）
nsteps                  = 5000        ; 最大步数（确保收敛）
emtol                   = 1000.0      ; 收敛判据（力的最大值≤1000 kJ·mol⁻¹·nm⁻¹）
emstep                  = 0.01        ; 步长（nm）

; 非键相互作用（与后续模拟一致）
cutoff-scheme           = Verlet
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0

; 输出设置（记录能量变化）
nstenergy               = 50          ; 每50步输出一次能量
energygrps              = System      ; 输出整个体系的能量
# 生成tpr文件
gmx grompp -f em.mdp -c 1s40_solv_ions.gro -p 1s40.top -o 1s40_em.tpr -maxwarn 1
# 运行能量最小化
gmx mdrun -v -deffnm 1s40_em
##平衡模拟
创建nvt.mdp
; 运行控制
integrator              = md   ; 温度耦合算法（稳定且高效）
nsteps                  = 100000      ; 总步数（100 ps，时间步长2 fs）
dt                      = 0.002       ; 时间步长（2 fs，需与约束匹配）

; 约束（固定蛋白/DNA重原子，加速平衡）
constraints             = h-bonds     ; 约束氢键（允许更大时间步长）
constraint_algorithm    = LINCS       ; 约束算法（适配蛋白/核酸）

; 温度耦合（300 K，弱耦合）
tcoupl                  = V-rescale
tc-grps                 = Protein_DNA  Water_and_ions  ; 分两组耦合（避免“热溶剂-冷溶质”）
tau_t                   = 0.1          0.1      ; 耦合时间常数（ps）
ref_t                   = 300          300      ; 参考温度（K）

; 非键相互作用
cutoff-scheme           = Verlet
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0

; 输出设置（每100 ps保存轨迹/能量）
nstxout                 = 50000       ; 每50000步（100 ps）保存结构
nstvout                 = 50000       ; 每50000步保存速度
nstenergy               = 5000        ; 每5000步（10 ps）保存能量
energygrps              = Protein_DNA  Water_and_ions
###创建index,group与nvt.mdp中的grps相同
gmx make_ndx -f 1s40_em.gro -o index.ndx
交互式窗口输入3｜1
q保存并退出
#####
#运行
gmx grompp -f nvt.mdp -c 1s40_em.gro -p 1s40.top -n index.ndx -o 1s40_nvt.tpr -maxwarn 1
gmx mdrun -v -deffnm 1s40_nvt
关键检查：用 gmx energy 查看温度曲线（输入 16 选择 Temperature），确认温度稳定在 300 K 附近。
gmx energy -f 1s40_nvt.edr -o 1s40_nvt.xvg
#创建npt.mdp
; 复用NVT参数，仅修改压力耦合部分
integrator              = md
nsteps                  = 200000      ; 200 ps（更长时间确保压力稳定）
dt                      = 0.002
constraints             = h-bonds
constraint_algorithm    = LINCS

; 温度耦合（与NVT一致）
tcoupl                  = V-rescale
tc-grps                 = Protein_DNA  Water_and_ions
tau_t                   = 0.1          0.1
ref_t                   = 300          300

; 压力耦合（1 bar，各向同性）
pcoupl                  = Parrinello-Rahman  ; 适合后续生产模拟的压力算法
pcoupltype              = isotropic          ; 各向同性压力耦合
tau_p                   = 2.0                ; 压力耦合时间常数（ps）
ref_p                   = 1.0                ; 参考压力（bar）
compressibility         = 4.5e-5             ; 水的压缩系数（bar⁻¹）

; 非键相互作用（与NVT一致）
cutoff-scheme           = Verlet
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0

; 输出设置
nstxout                 = 100000      ; 每100000步（200 ps）保存结构
nstvout                 = 100000
nstenergy               = 10000       ; 每10000步（20 ps）保存能量
energygrps              = Protein_DNA  Water_and_ions
gmx grompp -f npt.mdp -c 1s40_nvt.gro -t 1s40_nvt.cpt -p 1s40.top -o 1s40_npt.tpr -n index.ndx -maxwarn 1
gmx mdrun -v -deffnm 1s40_npt
##关键检查：用 gmx energy 查看压力曲线（输入 18 选择 Pressure）和密度曲线（输入 24 选择 Density），确认压力稳定在 1 bar、密度稳定在～1000 kg/m³（水的密度）。
gmx energy -f 1s40_npt.edr -o 1s40_npt.xvg
##生产模式
创建md.mdp
; 运行控制（以100 ns为例，nsteps = 100000000 / 2 = 50000000）
integrator              = md
nsteps                  = 50000000    ; 50000000步 × 2 fs/步 = 100 ns
dt                      = 0.002

; 约束与耦合（与NPT一致）
constraints             = h-bonds
constraint_algorithm    = LINCS
tcoupl                  = V-rescale
tc-grps                 = Protein_DNA  Water_and_ions
tau_t                   = 0.1          0.1
ref_t                   = 300          300
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5

; 非键相互作用（与NPT一致）
cutoff-scheme           = Verlet
rlist                   = 1.0
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0

; 输出设置（每100 ps保存轨迹，平衡后的数据更可靠）
nstxout                 = 50000       ; 50000步 × 2 fs = 100 ps，保存结构轨迹
nstvout                 = 0           ; 关闭速度输出（节省空间）
nstenergy               = 10000       ; 每20 ps保存能量
nstxtcout               = 50000       ; 每100 ps保存轨迹文件（.xtc）
xtc_precision           = 1000        ; 轨迹精度（1000 = 0.001 nm，平衡需求）
energygrps              = Protein_DNA  Water_and_ions
##生产模拟
gmx grompp -f md.mdp -c 1s40_npt.gro -t 1s40_npt.cpt -p 1s40.top -o 1s40_md.tpr -n index.ndx -maxwarn 1
gmx mdrun -v -deffnm 1s40_md -nt 8  ; -nt 25 用25核CPU运行（根据硬件调整）

可用 GROMACS 的分析工具（如 gmx rmsd、gmx hbond、gmx distance）对蛋白 - DNA 复合物的动态行为进行分析
# 计算 Protein_DNA 合并组的骨架原子 RMSD（先对齐再计算）
gmx rms -f 1s40_md.xtc -s 1s40_md.tpr -n index.ndx -o rmsd_total.xvg -fit rot+trans
交互界面操作：
第一步选择 对齐原子组（如 Backbone，蛋白质骨架 + DNA 骨架）；
第二步选择 计算原子组（同上，确保与对齐组一致，避免偏差）。
##计算蛋白质与 ssDNA 的最小距离
gmx mindist -f 1s40_md.xtc -s 1s40_md.tpr -n index.ndx -od min_dist.xvg -d 0.5
# 提取数据并计算结合持续时间占比（以0.5 nm为阈值）
awk '{if($1!~/#/ && $1!~/@/) print $1, $2}' min_dist.xvg > min_dist_data.txt
total=$(wc -l < min_dist_data.txt)
bound=$(awk '$2 < 0.5 {count++} END {print count}' min_dist_data.txt)
echo "结合持续时间占比：$((bound * 100 / total))%"

> 3  # 选择Protein组（第一组）
> 1  # 选择ssDNA组（第二组）

awk '/^Frame/ {frame=$1} $3 < 0.45 {print frame, $1, $2, $3}' all_dist.out > close_contacts.txt

下游处理
#####mindist,xvg --蛋白质与SSDNA的接触频率
touch mindist.sh
vi mindist.sh
#!/bin/bash

# ------------------------------
# 配置参数（仅需修改这里）
# ------------------------------
INPUT_XVG="min_dist.xvg"  # 输入的min_dist.xvg文件路径（相对或绝对路径）
OUTPUT_DIR="min_dist_results"  # 输出结果文件夹（自动创建）
THRESHOLD=0.5  # 结合阈值（nm），默认0.5 nm
PLOT_TITLE="蛋白质与ssDNA的最小距离随时间变化"  # 图表标题


# ------------------------------
# 自动处理流程
# ------------------------------
echo "===== 开始处理 $INPUT_XVG ====="

# 1. 创建输出文件夹
mkdir -p $OUTPUT_DIR
echo "输出结果将保存至：$(realpath $OUTPUT_DIR)"

# 2. 提取有效数据（去除注释行）
awk '{if($1!~/#/ && $1!~/@/) print $1, $2}' $INPUT_XVG > $OUTPUT_DIR/min_dist_clean.txt
if [ $? -ne 0 ]; then
    echo "错误：提取数据失败，请检查输入文件路径是否正确"
    exit 1
fi

# 3. 计算关键统计指标
total_frames=$(wc -l < $OUTPUT_DIR/min_dist_clean.txt)
if [ $total_frames -eq 0 ]; then
    echo "错误：未提取到有效数据，请检查输入文件格式"
    exit 1
fi

bound_frames=$(awk -v t=$THRESHOLD '$2 < t {count++} END {print count}' $OUTPUT_DIR/min_dist_clean.txt)
avg_dist=$(awk '{sum+=$2} END {printf "%.3f", sum/NR}' $OUTPUT_DIR/min_dist_clean.txt)
min_dist=$(awk 'NR==1{min=$2} $2<min{min=$2} END {printf "%.3f", min}' $OUTPUT_DIR/min_dist_clean.txt)
max_dist=$(awk 'NR==1{max=$2} $2>max{max=$2} END {printf "%.3f", max}' $OUTPUT_DIR/min_dist_clean.txt)
total_time=$(awk 'END {print $1/1000}' $OUTPUT_DIR/min_dist_clean.txt)  # 总时长（ns）
bound_ratio=$(echo "scale=2; $bound_frames * 100 / $total_frames" | bc)

# 4. 保存统计结果到文本文件
echo "===== 统计指标 =====" > $OUTPUT_DIR/statistics.txt
echo "总模拟时长：$total_time ns" >> $OUTPUT_DIR/statistics.txt
echo "总帧数：$total_frames" >> $OUTPUT_DIR/statistics.txt
echo "平均最小距离：$avg_dist nm" >> $OUTPUT_DIR/statistics.txt
echo "最小距离：$min_dist nm" >> $OUTPUT_DIR/statistics.txt
echo "最大距离：$max_dist nm" >> $OUTPUT_DIR/statistics.txt
echo "结合阈值：$THRESHOLD nm" >> $OUTPUT_DIR/statistics.txt
echo "结合持续时间占比：$bound_ratio%" >> $OUTPUT_DIR/statistics.txt
echo "统计结果已保存至：$OUTPUT_DIR/statistics.txt"

# 5. 调用Python脚本生成可视化图片
cat << EOF > $OUTPUT_DIR/plot_min_dist.py
import pandas as pd
import matplotlib.pyplot as plt
import os

# 读取数据
df = pd.read_csv(
    "$OUTPUT_DIR/min_dist_clean.txt",
    sep="\s+",
    header=None,
    names=["time_ps", "distance_nm"]
)
df["time_ns"] = df["time_ps"] / 1000  # 转换为ns

# 设置风格
plt.style.use("seaborn-v0_8-talk")
fig, ax = plt.subplots(figsize=(10, 5))

# 绘制最小距离曲线
ax.plot(
    df["time_ns"],
    df["distance_nm"],
    color="#2ecc71",
    linewidth=2,
    alpha=0.8,
    label="最小距离"
)

# 添加阈值线和平均线
ax.axhline(
    y=$THRESHOLD,
    color="#e74c3c",
    linestyle="--",
    linewidth=2,
    label=f"结合阈值 ({THRESHOLD} nm)"
)
ax.axhline(
    y=$avg_dist,
    color="#3498db",
    linestyle=":",
    linewidth=2,
    label=f"平均距离 ({avg_dist} nm)"
)

# 坐标轴与标题
ax.set_xlabel("时间 (ns)", fontsize=12)
ax.set_ylabel("最小距离 (nm)", fontsize=12)
ax.set_title("$PLOT_TITLE", fontsize=14, pad=20)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# 保存图片
output_img = os.path.join("$OUTPUT_DIR", "min_distance_plot.png")
plt.tight_layout()
plt.savefig(output_img, dpi=300, bbox_inches="tight")
print(f"可视化图片已保存至：{output_img}")
EOF

# 运行Python脚本
python3 $OUTPUT_DIR/plot_min_dist.py
if [ $? -ne 0 ]; then
    echo "警告：绘图失败，请检查是否安装pandas和matplotlib（可通过 pip install pandas matplotlib 安装）"
else
    echo "可视化图片生成成功：$OUTPUT_DIR/min_distance_plot.png"
fi

echo "===== 处理完成 ====="


####
chmod +x mindist.sh
./mindist.sh



#####pymol 动画轨迹可视化
gmx trjconv -s 1s40_md.tpr -f 1s40_md.xtc -o animated.pdb #不做中心去离子化，否则报错
pymol可视化


# 在PyMOL中加载处理后的轨迹文件
load "/Volumes/scyp1/细胞/cloudapi/animated.pdb"

# 设置显示样式
hide everything
show cartoon
color gray

# 可选：高亮特定区域（如配体）
select ligand, resn LIG  # 将LIG替换为实际的配体名称
show sticks, ligand
color red, ligand

# 创建旋转动画（60帧，每帧旋转6度）
mset 1 x60
mdo 1:60: rotate y, 6

# 播放预览
mplay