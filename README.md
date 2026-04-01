# ONLPoisson（Python 运行指南）

本仓库当前以 **Python 运行流程** 为主，用于组织 ONLPoisson 复现实验的执行入口与结果落盘。

> 说明：仓库中的 `.m` / `.c` / `.mexw64` / `.dll` 文件仅作历史参考，不参与当前复现实验执行链路。

## 1. 环境要求

- Python 3.10+
- 建议操作系统：Linux / macOS / Windows
- 依赖见 `requirements.txt`

## 2. 安装步骤

### 方式 A：使用 pip 安装（推荐）

```bash
python -m pip install -U pip
python -m pip install -r requirements.txt
python -m pip install -e .
```

### 方式 B：仅安装运行依赖

```bash
python -m pip install -U pip
python -m pip install -r requirements.txt
```

## 3. 运行命令

### 命令行入口（安装 `-e .` 后）

```bash
onlpoisson-run --images images.mat --output-dir results
```

### 或使用模块方式

```bash
python -m onlpoisson.cli --images images.mat --output-dir results
```

## 4. 结果文件位置

默认输出目录：`results/`

运行后会生成：

- `results/run_info.json`：本次运行的基础元数据（输入文件、图像键名、shape 等）

## 5. 历史文件说明

以下文件保留用于历史代码参考与论文实现对照：

- MATLAB 脚本与函数：`.m`
- C/MEX 源码与二进制：`.c` / `.mexw64` / `.dll`

这些文件 **不属于当前 Python 复现实验的执行链路**。
