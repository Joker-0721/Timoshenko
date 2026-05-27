import vtk
import numpy as np
from pathlib import Path
import sys

# 強制輸出編碼為 UTF-8（解決 Windows 終端機編碼問題）
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def read_vtu_modes(filepath, max_modes=50):
    """讀取 VTU 檔案中的所有 w1, w2, ... 陣列"""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(filepath))
    reader.Update()
    grid = reader.GetOutput()
    point_data = grid.GetPointData()
    
    modes = []
    for i in range(1, max_modes+1):
        arr = point_data.GetArray(f"w{i}")
        if arr is None:
            break
        modes.append(np.array(arr))
    return modes

def compute_mac_matrix(phi_ref, phi_comp):
    n_ref = len(phi_ref)
    n_comp = len(phi_comp)
    mac = np.zeros((n_ref, n_comp))
    for i in range(n_ref):
        for j in range(n_comp):
            a = phi_ref[i]
            b = phi_comp[j]
            num = (np.dot(a, b))**2
            den = np.dot(a, a) * np.dot(b, b)
            mac[i, j] = num / den if den != 0 else 0
    return mac

def main():
    exact_path = Path(r"D:\Joker\Timoshenko\vtk\mindlin_2d_ssss_q4_17x17_exact_modes.vtu")
    fem_path = Path(r"D:\Joker\Timoshenko\vtk\mindlin_2d_ssss_q4_17x17_vibration_FEM.vtu")
    
    if not exact_path.exists() or not fem_path.exists():
        print("ERROR: File not found. Check paths.")
        print("Exact path:", exact_path)
        print("FEM path:", fem_path)
        return
    
    print("Reading exact modes...")
    exact_modes = read_vtu_modes(exact_path, max_modes=50)
    print(f"Read {len(exact_modes)} exact modes")
    
    print("Reading FEM modes...")
    fem_modes = read_vtu_modes(fem_path, max_modes=50)
    print(f"Read {len(fem_modes)} FEM modes")
    
    # 取前 31 個比較
    n = min(31, len(exact_modes), len(fem_modes))
    exact_modes = exact_modes[:n]
    fem_modes = fem_modes[:n]
    
    print(f"\nComputing {n}x{n} MAC matrix ...")
    mac = compute_mac_matrix(exact_modes, fem_modes)
    
    # 儲存 CSV
    csv_path = "mac_17x17_31x31.csv"
    np.savetxt(csv_path, mac, delimiter=",", fmt="%.4f")
    print(f"MAC matrix saved to {csv_path}")
    
    # 顯示對角線
    print("\nDiagonal MAC values (exact vs FEM same index):")
    for i in range(n):
        print(f"exact w{i+1:2d} vs FEM w{i+1:2d}: {mac[i,i]:.4f}")
    
    # 最佳匹配
    print("\nBest match for each exact mode (MAC > 0.5):")
    for i in range(n):
        best_j = np.argmax(mac[i, :])
        best_val = mac[i, best_j]
        if best_val > 0.5:
            print(f"exact w{i+1:2d} -> FEM w{best_j+1:2d}  (MAC = {best_val:.4f})")
        else:
            print(f"exact w{i+1:2d} -> no clear match (best MAC = {best_val:.4f})")
    
    # 檢查成對區塊對角線和
    print("\n2x2 block diagonal sums (should be 1.0 for mode pairs):")
    for i in range(0, n-1, 2):
        if i+1 < n:
            block_sum = mac[i, i] + mac[i+1, i+1]
            print(f"Block w{i+1}-w{i+2} : {block_sum:.3f}")

if __name__ == "__main__":
    main()