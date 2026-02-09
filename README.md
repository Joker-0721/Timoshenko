# Timoshenko

本仓库整理与实现**鐵摩辛柯梁（Timoshenko beam）理论**在屈曲（buckling）问题中的若干关键主题，重点包括：

- 剪力变形与截面转角共同作用下的梁模型（相对 Euler–Bernoulli 梁更适用于厚梁/短粗梁等情况）
- 有限元素离散中的**剪力锁定（shear locking）**与“无锁定（locking-free）”思路
- 与屈曲临界载荷相关的线性特征值问题构造与数值求解
- Engesser 与 Haringx 两种常见修正/推导脉络的笔记整理

## 仓库内容

### 文档（Markdown）

- engesser.md：面向入门的“无锁定（locking-free）”概念梳理，包含剪力锁定的成因、影响，以及常见缓解方法（降阶积分、混合式、假设应变/应力等）与其在屈曲分析中的意义。
- timoshenko/：更偏推导与资料汇总的笔记目录，包含 Engesser/Haringx 相关文档与图示素材。

### 代码（Julia）

src/ 目录提供用于验证/演示的线性屈曲算例脚本，核心目标是构建并求解屈曲临界载荷对应的特征值问题。

- src/cantilever_buckling_linear.jl：悬臂梁（cantilever）线性屈曲示例。脚本中定义材料/几何参数（E、ν、L、b、h、I、A、κ、G 等），并以一定的单元数进行离散，组装相应矩阵后求解临界值。
- src/two_member_internal_hinge_fixed_fixed_buckling_linear.jl：两构件、内铰、两端固支（fixed-fixed）配置的线性屈曲示例。

## 关键概念（简述）

### 1) 为什么用鐵摩辛柯梁
鐵摩辛柯梁用两个基本未知量描述：横向位移 $w(x)$ 与截面转角 $\phi(x)$。它显式考虑剪切变形，因此当构件不够细长、或剪切效应不可忽略时，通常比 Euler–Bernoulli 梁更贴近真实响应。

### 2) 剪力锁定（Shear Locking）
在细长梁极限下，真实解趋向“纯弯曲”（剪切应变趋近于 0）。若使用标准低阶位移型单元直接离散鐵摩辛柯梁，离散空间可能无法同时满足该一致性约束，导致数值解表现得**过度刚硬**，位移被“锁住”，屈曲临界载荷也可能被显著高估。

### 3) 无锁定（Locking-Free）思路
常见策略包括：

- 降阶积分：对剪切相关项采用更低阶积分规则（需注意可能引入零能量模态等副作用）。
- 混合式/变分原理：将位移与应力/应变视为独立未知量（如基于 Hellinger–Reissner）。
- 假设应变/应力法：在单元内引入满足一致性条件的应变/应力场假设，提升薄梁极限下的表现。

## 如何运行（Julia）

前置条件：安装 Julia（建议 1.9+ 或更高版本）。脚本使用标准库 LinearAlgebra。

在仓库根目录运行：

```bash
julia --project=. src/cantilever_buckling_linear.jl
```

或：

```bash
julia --project=. src/two_member_internal_hinge_fixed_fixed_buckling_linear.jl
```

如果你没有使用 Julia project 环境，也可以直接运行（不依赖额外包时通常可行）：

```bash
julia src/cantilever_buckling_linear.jl
```

## 目录结构

```text
.
├── README.md
├── engesser.md
├── src/
│   ├── cantilever_buckling_linear.jl
│   └── two_member_internal_hinge_fixed_fixed_buckling_linear.jl
```

## 说明

- 本仓库当前更偏“推导笔记 + 算例脚本”，并非打包后的通用有限元库。
- 若你希望我把 README 再补上“每个脚本的输入/输出含义、关键矩阵与边界条件对应关系、以及典型结果如何解读”，告诉我你想以哪个脚本为主（cantilever 或 two-member）。