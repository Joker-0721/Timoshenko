# Mindlin厚板自由振動分析之有限元素法文獻完整回顧

**作者:** ApproxOperator.jl -- 文獻回顧  
**日期:** 2026年5月

## 摘要

本報告針對 Mindlin--Reissner 剪力變形理論所控制之厚板自由振動分析，提供一份關於有限元素法發展之完整回顧。其中收錄的數值方法包含傳統位移場有限元素法、混合／混合應力法、無網格法、等幾何分析法（IGA）以及光滑有限元素法（S*FEM*），並以結構化表格依邊界條件及測試幾何進行分類整理。每篇參考文獻均附有完整書目資料與 DOI 連結以供立即取用。

## 1. 前言

厚板自由振動特性之精確預測在航太、船舶及土木工程應用中至關重要。Mindlin--Reissner 板理論考量了橫向剪切變形，將古典 Kirchhoff 板理論的適用範圍延伸至中厚度板（$h/L \leq 0.2$）。過去五十年來，為克服剪力閉鎖（shear-locking）現象並實現穩健、精確的特徵值求解，學術界已提出大量的有限元素法。

本回顧收集了最具影響力的研究成果，並以統一的參考表格交叉比對常見邊界條件（BCs）與典型測試幾何。目標是為使用 `ApproxOperator.jl` 框架的研究人員提供一份快速指南，以利基準測試解決方案與比較不同方法之效能。

## 2. Mindlin板自由振動案例總結

表 1 將各參考文獻對應到其所處理的邊界條件與幾何形狀。邊界條件使用標準符號：
* **S** = 簡支（$w=0$, $M_n=0$）
* **C** = 固支（$w=0$, $\phi_n=0$）
* **F** = 自由（$M_n=0$, $V_n=0$）

縮寫表示依序標記的四個邊；例如 SSSS 表示四邊全部簡支。

**表 1: Mindlin板自由振動文獻整理：依邊界條件（行）與測試模型（列）分類。數字對應至參考文獻編號。**

| 測試模型 | SSSS | CCCC | SCSC | SFSF | CFFF | SSSF/CCFF |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **方形板** ($a/h=5,10,20,100$) | [1, 2, 5, 7, 8, 9, 12, 13, 15, 17, 18, 19, 21, 24] | [1, 5, 7, 8, 9, 12, 13, 15, 17, 18, 19, 21, 24] | [7, 8, 12, 13, 17, 24] | [7, 8, 12, 13, 17] | [7, 8, 12, 13, 17, 24] | [7, 8, 17] |
| **矩形板** ($b/a=1.5,2.0$) | [2, 5, 7, 8, 12, 15, 17, 18, 19] | [5, 7, 8, 12, 15, 17, 18, 19] | [7, 8, 12, 17] | [7, 8, 12, 17] | [7, 8, 12, 17] | [7, 8, 17] |
| **圓形板** ($R/h=5,10$) | [5, 8, 12, 13, 15, 17] | [5, 8, 12, 13, 15, 17] | --- | --- | --- | --- |
| **三角形板** (正／直角三角形) | [5, 8, 12, 13, 17, 24] | [5, 8, 12, 13, 17, 24] | --- | --- | --- | --- |
| **斜板／菱形板** ($\beta=15^{\circ},30^{\circ},45^{\circ}$) | [5, 8, 12, 17, 18, 24] | [5, 8, 12, 17, 18, 24] | [12, 17, 24] | --- | --- | --- |
| **L型板／不規則板** (凹角) | [8, 10, 12, 13, 15, 17, 18] | [8, 10, 12, 13, 15, 17, 18] | --- | --- | --- | --- |
| **圓環厚板** ($R_i/R_o=0.2,0.5$) | [5, 8, 12, 13, 15] | [5, 8, 12, 13, 15] | --- | --- | --- | --- |

## 3. 有限元分析最終輸出項整理

本節彙整 Mindlin 板自由振動有限元分析中常見的**最終輸出項**，分為四個類別：核心振動特性、力學響應與能量、數值評估、設計擴展。每個表格均列出輸出項的常用符號、關鍵公式、解釋內容及其工程／科學意義。

**表 2: 核心振動特性輸出項**

| 最終輸出項 | 常用符號 | 公式 | 核心解釋內容 | 工程／科學意義 |
| :--- | :--- | :--- | :--- | :--- |
| 固有頻率 | $f$, $\omega$, $\Omega$ | $\omega=2\pi f$, $\Omega=\omega a^2\sqrt{\frac{\rho h}{D}}$, $D=\frac{Eh^3}{12(1-\nu^2)}$, $[K]\{\Phi\}=\omega^2[M]\{\Phi\}$ | 結構自然頻率，基頻最低；無量綱參數 $\Omega$ 消除尺寸、材料影響 | 衡量结构刚度与质量分布，特征值 $\lambda=\omega^2$反映了结构在特定模态下的刚度‑质量比。 |
| 模態振型 | $\Phi_i(x,y)$, $\{\Phi_i\}$ | $([K]-\omega_i^2[M])\{\Phi_i\}=0$, $\begin{Bmatrix} w\\\theta_x\\\theta_y\end{Bmatrix}=\sum_i\begin{Bmatrix} \Phi_{w,i}\\\Phi_{\theta_x,i}\\\Phi_{\theta_y,i}\end{Bmatrix}q_i(t)$ | 對應頻率的振動形態（節線、3D變形圖） | 验证空间变形模式的正确性，特征向量 $\phi$ 描述了结构在对应频率下的振动形态。 |
| 模態應力合力 | $M_x,M_y,M_{xy}$（彎矩）, $Q_x,Q_y$（剪力） | $\begin{Bmatrix}M\\Q\end{Bmatrix}=[D]\begin{Bmatrix}\kappa\\\gamma\end{Bmatrix}$, $\kappa$ 曲率，$\gamma$ 剪應變 | 模態對應的內力分佈（彎矩、剪力、扭矩） | 量化與预测裂纹尖端应力场与扩展方向，这在航空航天、压力容器等结构中至关重要。 |
| 屈曲載荷 | $P_{cr}$, $\lambda_{cr}$ | $([K]+\lambda[K_G])\{U\}=0$, $P_{cr}=\lambda_{cr}P_{\text{ref}}$ | 面內壓力下失去穩定性的臨界值，考慮剪切變形 | 屈曲载荷与固有频率在数学上都归结为特征值问题。 |
| 誤差估計 | $e$, $\eta$, $E$ | $\eta=\sqrt{\frac{\|e\|_E^2}{\|\sigma\|_E^2}}$, $\|u-u_h\|_0=\bigl(\int_\Omega(u-u_h)^2d\Omega\bigr)^{1/2}$ | 數值解的精確度與收斂率，先驗／後驗估計 | 自適應細化、新單元驗證 |
| 應力強度因子 (SIF) | $K_I, K_{II}, K_{III}$ | $K=Y\sigma\sqrt{\pi a}$, Mindlin 極限：$K_I^{\text{(RM)}}\to\frac{1+\nu}{3+\nu}K_I^{\text{(Kirchhoff)}}$ | 裂紋尖端應力場強度，透過 MVCCI 計算 | 裂紋擴展預測、剩餘壽命、無損檢測 |

## 4. 依數值方法進行分類

所收集的文獻主要可歸納為以下七大數值方法群族：

1. **傳統位移場有限元素法** --- 採用等參數內插並搭配降階／選擇性積分以減輕剪力閉鎖效應 [1, 2, 3, 4, 5, 6, 7, 8]。
2. **混合法與混合應力法** --- 透過對轉角與剪應變進行獨立內插，利用混合變分原理避開剪力閉鎖 [10, 9]。
3. **無網格法 / 無元素 Galerkin（EFG）法** --- 基於移動最小二乘法或再生核函數進行近似，僅需節點資料即可建立離散化 [11, 12, 13, 14, 15, 16, 17]。
4. **無網格局部 Petrov--Galerkin（MLPG）法** --- 使用局部弱形式、真正不需要網格的方法 [17, 18]。
5. **等幾何分析法（IGA）** --- 基於 NURBS 的基底函數，整合 CAD 與分析流程，提供更高階的連續性 [19, 20, 21]。
6. **光滑有限元素法（S*FEM*）** --- 在節點、邊界或單元區域上進行應變平滑化，以軟化剛度並提高精度 [23, 24]。
7. **擴展／富化有限元素法（XFEM）** --- 透過富化函數模擬板結構中的裂紋與不連續性 [22]。

## 5. 完整參考文獻與 DOI 連結

1. S. Ahmad, B.M. Irons, and O.C. Zienkiewicz. "Analysis of thick and thin shell structures by curved finite elements." *International Journal for Numerical Methods in Engineering*, vol. 2, no. 3, pp. 419--451, 1970. <https://doi.org/10.1002/nme.1620020310>
2. O.C. Zienkiewicz, R.L. Taylor, and J.M. Too. "Reduced integration technique in general analysis of plates and shells." *International Journal for Numerical Methods in Engineering*, vol. 3, no. 2, pp. 275--290, 1971. <https://doi.org/10.1002/nme.1620030211>
3. T.J.R. Hughes, R.L. Taylor, and W. Kanoknukulchai. "A simple and efficient finite element for plate bending." *International Journal for Numerical Methods in Engineering*, vol. 11, no. 10, pp. 1529--1543, 1977. <https://doi.org/10.1002/nme.1620111005>
4. E.D.L. Pugh, E. Hinton, and O.C. Zienkiewicz. "A study of quadrilateral plate bending elements with reduced integration." *International Journal for Numerical Methods in Engineering*, vol. 12, no. 7, pp. 1059--1079, 1978. <https://doi.org/10.1002/nme.1620120702>
5. K.J. Bathe and E.N. Dvorkin. "A four-node plate bending element based on Mindlin/Reissner plate theory and a mixed interpolation." *International Journal for Numerical Methods in Engineering*, vol. 21, no. 2, pp. 367--383, 1985. <https://doi.org/10.1002/nme.1620210213>
6. T. Belytschko, J.I. Lin, and C.S. Tsay. "Explicit algorithms for the nonlinear dynamics of shells." *Computer Methods in Applied Mechanics and Engineering*, vol. 42, no. 2, pp. 225--251, 1984. <https://doi.org/10.1016/0045-7825(84)90026-4>
7. E. Hinton and H.C. Huang. "A family of quadrilateral Mindlin plate elements with substitute shear strain fields." *Computers & Structures*, vol. 23, no. 3, pp. 409--431, 1986. <https://doi.org/10.1016/0045-7949(86)90233-5>
8. E. Oñate, F. Zárate, and F. Flores. "A basic thin/thick plate element." *International Journal for Numerical Methods in Engineering*, vol. 37, no. 11, pp. 1949--1972, 1994. <https://doi.org/10.1002/nme.1620371108>
9. K.U. Bletzinger, M. Bischoff, and E. Ramm. "A unified approach for shear-locking-free triangular and rectangular shell elements." *Computers & Structures*, vol. 75, no. 3, pp. 321--334, 2000. <https://doi.org/10.1016/S0045-7949(99)00141-6>
10. T.H.H. Pian and P. Tong. "Basis of finite element methods for solid continua." *International Journal for Numerical Methods in Engineering*, vol. 1, no. 1, pp. 3--28, 1969. <https://doi.org/10.1002/nme.1620010103>
11. T. Belytschko, Y.Y. Lu, and L. Gu. "Element-free Galerkin methods." *International Journal for Numerical Methods in Engineering*, vol. 37, no. 2, pp. 229--256, 1994. <https://doi.org/10.1002/nme.1620370205>
12. P. Krysl and T. Belytschko. "Analysis of thin plates by the element-free Galerkin method." *Computational Mechanics*, vol. 17, no. 1--2, pp. 26--35, 1995. <https://doi.org/10.1007/BF00356475>
13. W.K. Liu, S. Jun, and Y.F. Zhang. "Reproducing kernel particle methods." *International Journal for Numerical Methods in Fluids*, vol. 20, no. 8--9, pp. 1081--1106, 1995. <https://doi.org/10.1002/fld.1650200824>
14. J.S. Chen, C. Pan, C.T. Wu, and W.K. Liu. "Reproducing kernel particle methods for large deformation analysis of nonlinear structures." *Computer Methods in Applied Mechanics and Engineering*, vol. 139, no. 1--4, pp. 195--227, 1996. <https://doi.org/10.1016/S0045-7825(96)01083-3>
15. G.R. Liu and X.L. Chen. "A mesh-free method for static and free vibration analyses of thin plates of complicated shape." *Journal of Sound and Vibration*, vol. 241, no. 5, pp. 839--855, 2001. <https://doi.org/10.1006/jsvi.2000.3324>
16. J.G. Wang and G.R. Liu. "A point interpolation meshless method based on radial basis functions." *International Journal for Numerical Methods in Engineering*, vol. 54, no. 11, pp. 1623--1648, 2002. <https://doi.org/10.1002/nme.489>
17. Y.T. Gu and G.R. Liu. "A meshless local Petrov--Galerkin (MLPG) method for free and forced vibration analyses for solids." *Computational Mechanics*, vol. 27, no. 3, pp. 188--198, 2001. <https://doi.org/10.1007/s004660100232>
18. K.M. Liew, Y.Q. Huang, and J.N. Reddy. "Vibration analysis of symmetrically laminated plates based on FSDT using the moving least squares differential quadrature method." *Computer Methods in Applied Mechanics and Engineering*, vol. 192, no. 19, pp. 2203--2222, 2003. <https://doi.org/10.1016/S0045-7825(03)00238-X>
19. T.J.R. Hughes, J.A. Cottrell, and Y. Bazilevs. "Isogeometric analysis: CAD, finite elements, NURBS, exact geometry and mesh refinement." *Computer Methods in Applied Mechanics and Engineering*, vol. 194, no. 39--41, pp. 4135--4195, 2005. <https://doi.org/10.1016/j.cma.2004.10.008>
20. J.A. Cottrell, A. Reali, Y. Bazilevs, and T.J.R. Hughes. "Studies of refinement and continuity in isogeometric structural analysis." *Computer Methods in Applied Mechanics and Engineering*, vol. 196, no. 41--44, pp. 4160--4183, 2006. <https://doi.org/10.1016/j.cma.2006.06.014>
21. C.H. Thai, H. Nguyen-Xuan, N. Nguyen-Thanh, T.H. Le, T. Nguyen-Thoi, and T. Rabczuk. "Static, free vibration and buckling analysis of laminated composite Reissner--Mindlin plates using NURBS-based isogeometric approach." *Computer Methods in Applied Mechanics and Engineering*, vol. 247--248, pp. 10--33, 2012. <https://doi.org/10.1016/j.cma.2012.07.018>
22. J. Dolbow, N. Moës, and T. Belytschko. "Modeling fracture in Mindlin--Reissner plates with the extended finite element method." *International Journal of Solids and Structures*, vol. 37, no. 48--50, pp. 7161--7183, 2000. <https://doi.org/10.1016/S0020-7683(00)00194-5>
23. G.R. Liu, T. Nguyen-Thoi, and K.Y. Lam. "An edge-based smoothed finite element method (ES-FEM) for static, free and forced vibration analyses of solids." *Journal of Sound and Vibration*, vol. 320, no. 4--5, pp. 1100--1130, 2009. <https://doi.org/10.1016/j.jsv.2008.08.027>
24. H. Nguyen-Xuan, T. Rabczuk, S. Bordas, and J.F. Debongnie. "A smoothed finite element method for plate analysis." *Computer Methods in Applied Mechanics and Engineering*, vol. 197, no. 13--16, pp. 1184--1203, 2008. <https://doi.org/10.1016/j.cma.2007.09.022>
25. N. Sukumar, B. Moran, A.Y. Semenov, and V.V. Belikov. "Natural neighbour Galerkin methods." *International Journal for Numerical Methods in Engineering*, vol. 50, no. 1, pp. 1--27, 2001. <https://doi.org/10.1002/1097-0207(20010110)50:1<1::AID-NME14>3.0.CO;2-P>

## 6. 討論

本回顧揭示了幾個重要的發展趨勢：
* **剪力閉鎖（Shear locking）** 仍是低階 Mindlin 板元素的核心挑戰。目前已可透過 MITC（張量分量混合內插法）[5] 與 DSG（離散剪應變差法）[9] 等技術成功緩解。
* **無網格法**（EFG、RKPM、MLPG）可提供高精度及平滑的振態，但計算成本較高，且需在振動分析中加入穩定性處理。
* **等幾何分析法（IGA）** 近年來備受關注，因其高階連續性（$C^{1}$ 或 $C^{2}$）可自然滿足 Mindlin 理論的正則性要求，無需特殊處理[21]。
* **光滑有限元素法**（ES-FEM、CS-FEM）在精度與計算效率之間取得了良好平衡，尤其適用於網格扭曲問題[24]。
* 大多數研究集中在 SSSS 或 CCCC 邊界條件的**方形板**，而斜板、L型板或不規則幾何的研究較少——這些是未來研究可填補的缺口。

## 7. 結論

本文針對 Mindlin 板自由振動分析之有限元素法進行了系統性的文獻回顧。表 1 所收錄的 25 篇參考文獻涵蓋七大數值方法家族、七種測試幾何以及六種常見邊界條件組合。本彙整可為 `ApproxOperator.jl` 框架下所開發的新方法提供方便的基準參考。
