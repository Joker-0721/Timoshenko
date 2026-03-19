# Lean4 安裝指南（Windows 10/11）

## 前置要求

在安裝 Lean4 之前，您需要先安裝以下軟體：

### 1. 安裝 Git
如果您還沒有 Git，請從以下網址下載並安裝：
- https://git-scm.com/download/win

安裝時建議使用預設選項。

### 2. 安裝 .NET SDK
Lean4 需要 .NET 環境。請從以下網址下載並安裝：
- https://dotnet.microsoft.com/download

選擇最新的 LTS 版本（.NET 8 或更新版本）。

---

## 安裝 Lean4

### 方法一：使用 elan（推薦）

elan 是 Lean 的版本管理器，安裝過程會自動處理所有依賴。

**步驟 1：打開 PowerShell**

按 `Win + X`，選擇「Windows PowerShell (系統管理員)」

**步驟 2：安裝 elan**

```powershell
curl --proto '=https' --tlsv1.2 -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh
```

**步驟 3：重啟終端機**

關閉 PowerShell 並重新打開，讓環境變數生效。

**步驟 4：驗證安裝**

```powershell
lean --version
elan --version
```

應該會看到類似以下的輸出：
```
lean version 4.0.0
elan 1.4.5
```

---

## VS Code 擴展安裝

### 1. 安裝 VS Code
如果您還沒有 VS Code，請從以下網址下載並安裝：
- https://code.visualstudio.com/

### 2. 安裝 Lean4 擴展

1. 打開 VS Code
2. 按 `Ctrl + Shift + X` 打開擴展市場
3. 搜尋 `Lean4`（由 Lean Prover 團隊開發）
4. 點擊「安裝」按鈕

安裝完成後，您應該會看到 Lean4 的工具列出現在 VS Code 底部。

---

## 建立第一個 Lean4 專案

### 步驟 1：建立專案資料夾

打開終端機，切換到您的工作目錄：

```powershell
cd d:\Joker\Timoshenko
mkdir Lean4Timoshenko
cd Lean4Timoshenko
```

### 步驟 2：初始化 Lean4 專案

```powershell
lean --init
```

這會建立以下檔案結構：
```
Lean4Timoshenko/
├── LeanToolchain
├── lake-manifest.json
├── lakefile.lean
└── Main.lean
```

### 步驟 3：在 VS Code 中打開專案

```powershell
code .
```

或手動在 VS Code 中選擇「File → Open Folder」並選擇 `Lean4Timoshenko` 資料夾。

---

## 驗證 Lean4 環境

打開 `Main.lean` 檔案，應該會看到預設的 "Hello, world!" 範例。

在 VS Code 中：
- 按 `Ctrl + Enter` 可以執行/ evaluate 游標所在位置的程式碼
- 底部狀態列會顯示 Lean 伺服器狀態（應該是 ✅）

---

## 常見問題

### Q: 點擊程式碼沒有反應
A: 確認 VS Code 底部狀態列的 Lean 伺服器是否正常運行。如果顯示紅色錯誤，嘗試重新載入視窗（按 `Ctrl + Shift + P`，輸入「Reload Window」）。

### Q: 如何更新 Lean 版本？
A: 使用以下命令：
```powershell
lean --update
```

### Q: 出現 "elan not found" 錯誤
A: 請重新啟動終端機，讓環境變數生效。

---

## 下一步

現在您可以開始學習 Lean4 基礎語法，並使用我們的 Timoshenko 梁理論教學例子！

請參考：[LEAN4_TUTORIAL.md](./LEAN4_TUTORIAL.md)
