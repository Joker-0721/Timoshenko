# ==============================================================================
#  高級計算固體力學研究平台：Mindlin 屈曲分析理論解析解對照組 (Exact)
#  優化：1. 完全移除 L2 誤差計算
#        2. 大循環 (15➔25) 數據實時追加至全局統一 CSV 大表 (method="Exact", 誤差恆為0)
# ==============================================================================

const LOCAL_APPROX_OPERATOR = normpath(joinpath(@__DIR__, "..", "..", "ApproxOperator.jl"))
if isdir(LOCAL_APPROX_OPERATOR) && !(LOCAL_APPROX_OPERATOR in LOAD_PATH)
    pushfirst!(LOAD_PATH, LOCAL_APPROX_OPERATOR)
end

using WriteVTK
using Printf

# 全局數據庫路徑規範
const DATA_DIR = "./data"
const UNIFIED_CSV = joinpath(DATA_DIR, "buckling_exact.csv")
mkpath(DATA_DIR)

println("="^80)
println(" 執行 Exact 模組：將理論基基準追加至全局大 CSV 表 (ndiv = 9 ➔ 25) ")
println("="^80)

for n_div in 9:25
    h_size = 1.0 / n_div
    log10_h = log10(h_size)
    
    # 理論經典簡支薄板解析解常數
    k_analytical = 4.0
    rel_error = 0.0
    dummy_lambda = 0.0  # 連續理論解無離散特徵值，記為 0.0
    
    # 🚀 數據流管線：追加寫入大 CSV 文件，確保 Exact 的基準也同步在內
    file_exists = isfile(UNIFIED_CSV)
    open(UNIFIED_CSV, "a") do io
        if !file_exists
            println(io, "method,ndiv,h,log10_h,lambda_cr,k_num,relative_error_k")
            file_exists = true
        end
        @printf(io, "Exact,%d,%.6e,%.6f,%.6e,%.6f,%.6e\n", 
                n_div, h_size, log10_h, dummy_lambda, k_analytical, rel_error)
    end
    println("  [Exact] ndiv = $(n_div) 理論基準已追加寫入大表。")
end