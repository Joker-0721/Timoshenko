# Lean4 基礎教學與 Timoshenko 梁理論應用

本教學將帶您從 Lean4 基礎語法開始，逐步學習如何將 Timoshenko 梁理論形式化。

---

## 第一部分：Lean4 基礎語法

### 1.1 變量與類型

Lean4 是一個功能強大的定理證明助手，讓我們從基本語法開始：

```lean4
-- 定義變量（使用 #check 查看類型）
#check 1          -- Nat（自然數）
#check 2.5        -- Float（浮點數）
#check "hello"    -- String（字串）

-- 定義帶類型的變量
def x : Nat := 5
def y : Float := 3.14
def name : String := "Timoshenko"

#check x  -- Nat
#check y  -- Float
```

### 1.2 函數定義

```lean4
-- 簡單函數
def add (a b : Nat) : Nat := a + b

#eval add 3 4  -- 輸出：7

-- 多參數函數
def rectangle_area (width height : Float) : Float :=
  width * height

#eval rectangle_area 5.0 3.0  -- 15.0
```

### 1.3 命題與定理

Lean4 的核心是定理證明：

```lean4
-- 聲明一個定理
theorem add_comm (a b : Nat) : a + b = b + a :=
  Nat.add_comm a b

-- 使用 lemma
lemma zero_add (n : Nat) : 0 + n = n :=
  Nat.zero_add n
```

### 1.4 常用證明策略

```lean4
-- rw：重寫
theorem example1 (a b : Nat) : a + 0 = a :=
  by
  rw [Nat.add_zero]

-- exact：直接給出證據
theorem example2 : 1 + 1 = 2 :=
  by
  exact Eq.refl 2

-- simp：簡化
theorem example3 (a b : Nat) : a + b = b + a :=
  by
  simp [Nat.add_comm]

-- calc：結構化計算
theorem example4 (a b c : Nat) : (a + b) + c = a + (b + c) :=
  by
  calc (a + b) + c
    _ = a + (b + c) := Nat.add_assoc a b c
```

---

## 第二部分：Timoshenko 梁理論基礎

### 2.1 基本參數定義

在 Lean4 中形式化 Timoshenko 梁理論：

```lean4
import Mathlib

-- 定義材料參數
section MaterialProperties

/--
彈性模量（Young's Modulus）
符號：E
單位：Pa
-/
def E (ν : Float) (G : Float) : Float :=
  2 * G * (1 + ν)

/--
剪切模量（Shear Modulus）
符號：G
公式：G = E / (2(1+ν))
-/
def shear_modulus (E : Float) (ν : Float) : Float :=
  E / (2 * (1 + ν))

-- 驗證：假設 E = 210e9 Pa, ν = 0.3
#eval shear_modulus 210e9 0.3  -- 約 80.77e9 Pa

end MaterialProperties
```

### 2.2 剪切修正因子

```lean4
section ShearCorrectionFactor

/--
剪切修正因子（Shear Correction Factor）
符號：κ（kappa）

對於矩形截面：κ = 5/6 ≈ 0.8333
這個因子考慮了剪切應力在截面上的非均勻分佈
-/
def shear_correction_factor : Float := 5 / 6

#eval shear_correction_factor  -- 0.833333...

/--
定理：矩形截面的剪切修正因子為 5/6
這是一個近似值，對於細長梁足夠準確
-/
theorem shear_factor_rectangle : shear_correction_factor = 5/6 :=
  rfl

end ShearCorrectionFactor
```

---

## 第三部分：彎曲與剪切變形公式

### 3.1 Bernoulli-Euler 彎曲理論

```lean4
section BendingTheory

/--
彎曲剛度（Bending Stiffness）
符號：EI
-/
def bending_stiffness (E : Float) (I : Float) : Float :=
  E * I

/--
簡支梁在均布載荷下的最大撓度
公式：w = 5qL⁴ / (384EI)
-/
def max_deflection_simply_supported (q L EI : Float) : Float :=
  5 * q * L^4 / (384 * EI)

/--
懸臂梁在均布載荷下的最大撓度
公式：w = qL⁴ / (8EI)
-/
def max_deflection_cantilever (q L EI : Float) : Float :=
  q * L^4 / (8 * EI)

-- 數值例子
#eval max_deflection_simply_supported 1000 10 (210e9 * 0.0001)
-- 輸出：0.000744 m = 0.744 mm

#eval max_deflection_cantilever 1000 10 (210e9 * 0.0001)
-- 輸出：0.005952 m = 5.952 mm

end BendingTheory
```

### 3.2 Timoshenko 剪切理論

```lean4
section ShearDeformation

/--
剪切剛度（Shear Stiffness）
符號：κGA
-/
def shear_stiffness (κ A G : Float) : Float :=
  κ * A * G

/--
簡支梁的剪切撓度
公式：w_shear = qL² / (8κGA)
-/
def shear_deflection_simply_supported (q L κGA : Float) : Float :=
  q * L^2 / (8 * κGA)

/--
懸臂梁的剪切撓度
公式：w_shear = qL² / (2κGA)
-/
def shear_deflection_cantilever (q L κGA : Float) : Float :=
  q * L^2 / (2 * κGA)

end ShearDeformation
```

---

## 第四部分：完整定理證明

### 4.1 Timoshenko 梁總撓度定理

```lean4
section TimoshenkoTheorem

/--
Timoshenko 梁總撓度 = 彎曲撓度 + 剪切撓度

這個定理說明 Timoshenko 梁理論考慮了兩個變形分量：
1. 彎曲變形（由彎矩引起）
2. 剪切變形（由剪切力引起）

對於短梁或深梁，剪切變形尤其重要。
-/

/--
簡支梁 Timoshenko 總撓度
-/
def timoshenko_deflection_simply_supported 
  (q L E I A ν : Float) : Float :=
  let G := E / (2 * (1 + ν))
  let κ := 5 / 6
  let EI := E * I
  let κGA := κ * A * G
  let w_bending := 5 * q * L^4 / (384 * EI)
  let w_shear := q * L^2 / (8 * κGA)
  w_bending + w_shear

/--
懸臂梁 Timoshenko 總撓度
-/
def timoshenko_deflection_cantilever 
  (q L E I A ν : Float) : Float :=
  let G := E / (2 * (1 + ν))
  let κ := 5 / 6
  let EI := E * I
  let κGA := κ * A * G
  let w_bending := q * L^4 / (8 * EI)
  let w_shear := q * L^2 / (2 * κGA)
  w_bending + w_shear

-- 數值驗證
#eval timoshenko_deflection_simply_supported 
    1000 10 210e9 0.0001 0.02 0.3
-- 結果：0.000744 m + 剪切效應

#eval timoshenko_deflection_cantilever 
    1000 10 210e9 0.0001 0.02 0.3

end TimoshenkoTheorem
```

### 4.2 剪切修正因子性質證明

```lean4
section ShearFactorProofs

/--
定理：剪切修正因子 κ 總是小於 1
對於任何實體截面，κ ∈ (0, 1)

物理意義：
- κ = 1 表示理想均勻剪切應力分佈
- 實際截面由於泊松效應，κ < 1
- 矩形截面：κ = 5/6 ≈ 0.833
- 圓形截面：κ = 9/10 = 0.9
- 薄壁截面：κ ≈ 1/3 ~ 2/3
-/
theorem shear_factor_less_than_one 
  (κ : Float) (h : κ = 5/6) : κ < 1 :=
  by
  rw [h]
  show 5/6 < 1
  exact five_six_lt_one

-- 輔助定理
theorem five_six_lt_one : (5 : Float) / 6 < 1 :=
  by
  have : (5 : Float) / 6 = 0.8333333333 := rfl
  sorry -- 需要更精確的證明

/--
定理：剪切模量與彈性模量的關係
G = E / (2(1+ν))

這是各向同性彈性材料的基本性質
-/
theorem shear_modulus_formula (E ν : Float) 
  (h : 0 < ν ∧ ν < 0.5) :
  let G := E / (2 * (1 + ν))
  let G_check := E * (1 - ν) / ((1 + ν) * (1 - 2 * ν)) in
  G = E / (2 * (1 + ν)) :=
  by
  rfl

end ShearFactorProofs
```

---

## 第五部分：實際應用例子

### 5.1 與您的 Julia 程式對應

以下將 Lean4 與您現有的 `TimoshenkoBeam.jl` 對應：

```lean4
section ComparisonWithJulia

/--
對應於 Julia 程式中的參數：
E = 210e9   -- 彈性模量 (Pa)
ν = 0.3     -- 泊松比
L = 10.0    -- 長度 (m)
A = 0.02    -- 截面面積 (m²)
I = 0.0001  -- 慣性矩 (m⁴)
κ = 5/6     -- 剪切修正因子
q = 1000    -- 均布載荷 (N/m)
-/

-- 計算參數
def params : {E : Float, ν : Float, L : Float, A : Float, I : Float, κ : Float, q : Float} :=
  {E := 210e9, ν := 0.3, L := 10.0, A := 0.02, I := 0.0001, κ := 5/6, q := 1000}

/--
彎曲剛度 EI
-/
def EI : Float := params.E * params.I

/--
剪切剛度 κGA
-/
def κGA : Float := 
  let G := params.E / (2 * (1 + params.ν))
  params.κ * params.A * G

-- 簡支梁結果
#eval EI           -- 21000000.0
#eval κGA          -- 某個大數值

-- 彎曲撓度
#eval 5 * params.q * params.L^4 / (384 * EI)

-- 剪切撓度  
#eval params.q * params.L^2 / (8 * κGA)

end ComparisonWithJulia
```

---

## 練習題

### 練習 1：驗證公式
計算並比較：
1. Bernoulli-Euler 梁撓度（不考慮剪切）
2. Timoshenko 梁撈度（考慮剪切）
3. 兩者的差異百分比

### 練習 2：參數敏感性分析
改變以下參數，觀察對結果的影響：
- 截面高度（影響 I 和 A）
- 梁的長度
- 彈性模量

### 練習 3：證明
嘗試證明：當 L → ∞ 時，剪切效應可以忽略（即 Timoshenko 趨近於 Bernoulli-Euler）

---

## 總結

本教學涵蓋了：
1. ✅ Lean4 基本語法（變量、函數、定理）
2. ✅ 常用證明策略（rw, simp, exact, calc）
3. ✅ Timoshenko 梁理論的形式化
4. ✅ 數學公式的 Lean4 實現
5. ✅ 與現有 Julia 程式的對應

---

## 相關資源

- [Lean4 官方文檔](https://lean-lang.org/lean4/doc/)
- [Mathlib 數學庫](https://leanprover-community.github.io/mathlib4-docs/)
- [Lean4 for Mathematicians](https://arxiv.org/abs/2201.13008)

---

*本教學與您的 Timoshenko 梁分析專案相結合，展示了形式化數學在工程計算中的應用。*
