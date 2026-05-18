# SSSS Mindlin 方板自由震動解析解對比

## 來源

- Xing, Y. and Liu, B. (2009), *Closed form solutions for free vibrations of rectangular Mindlin plates*, Acta Mechanica Sinica, 25, 689-698. https://pubs-en.cstam.org.cn/article/doi/10.1007/s10409-009-0253-7?viewType=HTML
- Huang and Huang (2020), Appendix A gives the simply supported Mindlin plate trigonometric series form used for SSSS separation. https://pmc.ncbi.nlm.nih.gov/articles/PMC7503694/

這裡使用的是 SSSS Mindlin/Reissner 厚板的 Navier/separation 型式。它和 FEM 模型一致地保留剪切變形與轉動慣量，所以比 Kirchhoff 薄板極限更適合 `h/b = 0.1` 的目前算例。

## 公式

對第 `(m,n)` 個解析候選模態：

```text
w = W sin(alpha x) sin(beta y)
theta_x = X cos(alpha x) sin(beta y)
theta_y = Y sin(alpha x) cos(beta y)

alpha = m*pi/a
beta  = n*pi/b
D = E*h^3/[12*(1-nu^2)]
G = E/[2*(1+nu)]
S = kappa*G*h
Ir = rho*h^3/12
```

解 3x3 廣義特徵值問題：

```text
(K_mn - omega^2 M_mn) q = 0
q = [W, X, Y]^T
M_mn = diag(rho*h, Ir, Ir)

K11 = S*(alpha^2 + beta^2)
K12 = -S*alpha
K13 = -S*beta
K22 = S + D*(alpha^2 + (1-nu)/2*beta^2)
K33 = S + D*(beta^2 + (1-nu)/2*alpha^2)
K23 = D*(1+nu)/2*alpha*beta
```

每組 `(m,n)` 有三個分支；本圖取最低 bending branch，排序後和 FEM 的 `mode_rank` 對比。

無因次頻率使用：

```text
Omega = omega*a^2*sqrt(rho*h/D)
```

## 為什麼要這樣

原圖的 residual 只能檢查 FEM 是否滿足離散後的 `K phi = omega^2 M phi`，不能判斷它是否接近連續方板解析解。加入解析頻率後，可以直接看出頻率偏硬或偏軟，也能檢查方板退化模態，例如 `(1,2)` 與 `(2,1)` 應有相同解析頻率。

## 目前參數

```text
E = 1e+08
nu = 0.3
a = 1
b = 1
h = 0.1
rho = 1
kappa = 0.83333333
D = 9157.5092
S = 3205128.2
Ir = 8.3333333e-05
```

退化檢查：`omega_exact(1,2) - omega_exact(2,1) = 0.000e+00`，符合方板解析解對稱性。

## 前 9 個模態對比

| mode | exact (m,n) | omega FEM | omega exact | freq FEM Hz | freq exact Hz | error % | boundary max | residual |
|---:|:---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | (1,1) | 10031.5 | 5769.32 | 1596.57 | 918.216 | +73.877 | 3.097e-10 | 7.214e-06 |
| 2 | (1,2) | 24816.7 | 13763.7 | 3949.7 | 2190.56 | +80.305 | 9.090e-10 | 1.537e-06 |
| 3 | (2,1) | 27693.6 | 13763.7 | 4407.57 | 2190.56 | +101.208 | 1.012e-09 | 1.045e-06 |
| 4 | (2,2) | 38741.2 | 21120.7 | 6165.85 | 3361.47 | +83.427 | 1.213e-09 | 8.848e-07 |
| 5 | (3,1) | 56531.5 | 25733.7 | 8997.27 | 4095.64 | +119.679 | 1.607e-09 | 1.615e-07 |
| 6 | (1,3) | 56994.7 | 25733.7 | 9070.98 | 4095.64 | +121.479 | 9.743e-10 | 3.113e-07 |
| 7 | (2,3) | 61588.1 | 32283.9 | 9802.05 | 5138.15 | +90.770 | 2.059e-09 | 3.618e-07 |
| 8 | (3,2) | 66316.3 | 32283.9 | 10554.6 | 5138.15 | +105.416 | 1.540e-09 | 1.411e-07 |
| 9 | (1,4) | 80917.2 | 40435.6 | 12878.4 | 6435.53 | +100.114 | 1.591e-09 | 1.473e-07 |

## 輸出檔

- `frequency_residual_boundary_with_mindlin_exact.png`
- `dimensionless_frequency_mindlin_exact_comparison.png`
- `2d_ssss_vibration_mindlin_exact_comparison.csv`
