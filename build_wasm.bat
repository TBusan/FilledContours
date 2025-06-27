@echo off
rem Script to build FilledContours as WebAssembly on Windows

echo ===== FilledContours WebAssembly构建脚本 =====
echo.

rem 设置Python路径
set PATH=D:\miniconda3;D:\miniconda3\Scripts;%PATH%

rem 设置MinGW路径
set PATH=D:\Programs\mingw64\bin;%PATH%
echo 当前PATH: %PATH%

rem 明确设置Emscripten环境
set EMSDK=D:\emsdk
set EM_CONFIG=%EMSDK%\.emscripten
set EMSCRIPTEN=%EMSDK%\upstream\emscripten
set PATH=%EMSCRIPTEN%;%EMSDK%;%EMSDK%\node\14.18.2_64bit\bin;%EMSDK%\upstream\bin;%PATH%

echo EMSDK = %EMSDK%
echo EMSCRIPTEN = %EMSCRIPTEN%
echo EM_CONFIG = %EM_CONFIG%
echo.

rem 确保已激活Emscripten环境
call %EMSDK%\emsdk_env.bat

rem 验证Emscripten环境
echo 验证Emscripten环境...
where emcc
if %ERRORLEVEL% NEQ 0 (
    echo [错误] emcc不在PATH中，无法继续。
    exit /b 1
)

emcc --version
echo.

rem 确保能找到mingw32-make
echo 验证MinGW环境...
where mingw32-make
if %ERRORLEVEL% NEQ 0 (
    echo [错误] mingw32-make不在PATH中，无法继续。
    exit /b 1
)
echo.

rem 删除并重新创建构建目录
echo 准备构建目录...
if exist "build-wasm" rmdir /s /q build-wasm
mkdir build-wasm
cd build-wasm
echo.

rem 配置CMake项目
echo 配置CMake项目...
echo emcmake cmake .. -G "MinGW Makefiles" -DCMAKE_MAKE_PROGRAM=mingw32-make -DCMAKE_BUILD_TYPE=Release -DCMAKE_CROSSCOMPILING_EMULATOR=node -DCMAKE_INSTALL_PREFIX=./install -DCMAKE_VERBOSE_MAKEFILE=ON -C ../CMakeLists_wasm.txt
call emcmake cmake .. -G "MinGW Makefiles" ^
    -DCMAKE_MAKE_PROGRAM=mingw32-make ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DCMAKE_CROSSCOMPILING_EMULATOR=node ^
    -DCMAKE_INSTALL_PREFIX=./install ^
    -DCMAKE_VERBOSE_MAKEFILE=ON ^
    -C ../CMakeLists_wasm.txt ^
    -DEMSCRIPTEN=1

if %ERRORLEVEL% NEQ 0 (
    echo [错误] CMake配置失败
    exit /b 1
)
echo.

rem 编译项目
echo 编译FilledContours WebAssembly...
call emmake mingw32-make -j4

if %ERRORLEVEL% NEQ 0 (
    echo [错误] 构建失败
    exit /b 1
)

echo.
echo ===== 构建成功! =====
echo 输出文件在 build-wasm\bin\ 目录中
echo.
echo 启动本地服务器:
echo cd build-wasm\bin ^&^& python -m http.server 8000
echo 然后访问浏览器: http://localhost:8000/
echo.

cd .. 