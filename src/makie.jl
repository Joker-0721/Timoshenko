using GLMakie, MeshIO, GeometryBasics
vtk = load("fem2.vtu")
mesh = Mesh(vtk)  # 自动转换为 Triangle 网格

data = vtk["v₁"]  # 直接取数组

fig, ax, plt = mesh(mesh, color=data, colormap=:jet, shading=false)
ax.title = "Eigenvector 1"
Colorbar(fig[1,2], plt)
fig