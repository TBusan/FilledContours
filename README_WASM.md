# FilledContours WebAssembly 构建指南

本文档提供了将FilledContours编译为WebAssembly并在浏览器中运行的详细说明。

## 先决条件

在开始之前，请确保已安装以下工具：

1. **Emscripten SDK (emsdk)**
   - 安装路径: `D:\emsdk`
   - 确保已激活最新版本: `emsdk install latest && emsdk activate latest`

2. **MinGW**
   - 安装路径: `D:\Programs\mingw64`
   - 确保mingw32-make在PATH中

3. **CMake**
   - 版本 3.15 或更高

4. **Python**
   - 用于运行本地HTTP服务器

5. **VTK (可选，用于完整版本)**
   - 需要使用Emscripten编译的VTK库
   - 默认路径: `D:\VTK-emscripten`

## 构建选项

本项目提供了几种不同的构建脚本，根据您的需求选择合适的脚本：

### 1. 简化版本 (无VTK依赖)

如果您只想测试基本的WebAssembly功能，可以使用简化版本：

```bash
.\build_wasm_simple.bat
```

这将创建一个简单的"Hello World"应用程序，输出文件位于`build-simple`目录中。

### 2. 完整版本 (需要VTK)

如果您想构建完整的FilledContours应用程序，包括VTK可视化功能：

```bash
.\build_wasm_final.bat
```

这将创建完整的FilledContours WebAssembly应用程序，输出文件位于`build-wasm\bin`目录中。

**注意**: 此选项需要预先使用Emscripten编译VTK库。如果您尚未编译VTK，请参考下面的说明。

### 3. 诊断工具

如果您在构建过程中遇到问题，可以使用以下诊断工具：

```bash
.\simplest_test.bat     # 测试基本的Emscripten + CMake环境
.\test_emscripten.bat   # 测试Emscripten环境
.\debug_build.bat       # 详细的构建诊断
```

## 运行WebAssembly应用程序

构建完成后，您可以使用以下命令启动本地HTTP服务器：

```bash
.\run_wasm.bat
```

然后在浏览器中访问：
- 简化版本: `http://localhost:8000/FilledContours_simple.html`
- 完整版本: `http://localhost:8000/`

## 编译VTK为WebAssembly (可选)

如果您需要构建完整版本，但尚未编译VTK库，请按照以下步骤操作：

1. 克隆VTK源代码：
```bash
git clone https://github.com/Kitware/VTK.git vtk-src
cd vtk-src
git checkout v9.2.0  # 或其他稳定版本
```

2. 使用提供的脚本编译VTK：
```bash
.\setup_vtk_emscripten.bat
```

这将在`D:\VTK-emscripten`目录中创建编译好的VTK库。

## 故障排除

如果您在构建过程中遇到问题，请检查以下事项：

1. 确保Emscripten已正确安装并激活
2. 确保MinGW在PATH中
3. 如果使用完整版本，确保VTK已使用Emscripten编译
4. 检查日志输出以获取详细错误信息

## 文件结构

- `build_wasm_simple.bat` - 构建简化版本
- `build_wasm_final.bat` - 构建完整版本
- `run_wasm.bat` - 启动HTTP服务器
- `simplest_test.bat` - 测试基本环境
- `debug_build.bat` - 详细诊断
- `setup_vtk_emscripten.bat` - 编译VTK库
- `FilledContours1.cxx` - 主源代码
- `shell.html` - HTML模板
- `web/` - Web界面文件
  - `index.html` - 主页面
  - `style.css` - 样式表
  - `filledcontours-api.js` - JavaScript API

## 注意事项

- WebAssembly应用程序必须通过HTTP服务器访问，直接打开HTML文件将不起作用
- 确保浏览器支持WebAssembly和WebGL2
- 大型数据集可能需要增加内存限制，可以在编译标志中调整 