# ONLPoisson 最小改动复现计划

## 0. 目标和边界

目标：在原 repo 基础上，用**最小改动**复现这篇 paper 自己的方法（也就是 Table 1 里的 **Ours**）的实验流程与数值趋势。

边界：

1. 第一阶段**只复现 Ours**，不在这个 repo 里顺手补全所有 baseline。
2. 换语言，把 MATLAB/C 整体改写成 Python。
3. 不做大重构，不做“顺手优化”。
4. 能复用原 repo 的地方就复用；只有和论文明显不一致的地方才改。

一句话原则：**保留 repo 的骨架，只修正和论文不对齐的部分。**

---

## 1. 实现语言与依赖边界（统一为 Python）

1. 目标实现语言统一为 **Python**。
2. MATLAB/C 文件（如 `fnmise.m`、`nlm_poisson004.c`、MEX 相关产物）仅作为公式与实现细节参考，不再作为执行依赖。
3. 复现实验脚本统一使用 Python，覆盖完整流程：数据加载、噪声生成、去噪与统计。
4. 明确输出物为：Python 源码、可运行脚本、结果表格（mean/std NMISE）。

说明：参数与实验协议保持不变（图像集、peaks、30 次试验、固定 Gaussian Step 2）。

---

## 2. 哪些地方必须改

## 2.1 不再用原 `demo_onl.m` 直接复现论文

最小方案：**新增一个新的 runner**，例如：

- `reproduce_table1_ours.m`
- `build_mex.m`
- （可选）`load_paper_images.m`

不要在旧 `demo_onl.m` 上继续堆大量 `if/else`。新增 runner 的改动更小、也更干净。

### 原 `demo_onl.m` 的问题

1. 它顶部注释还在说旧的 LNLA2009 代码，不是 2021 这篇 paper 的实验脚本。
2. 它现在用的图像集合是 `Spots, Galaxy, Ridges, Barbara, Cells`，和 paper 的 `Spots, Galaxy, Cells, Texture, Fermi` 不一致。
3. 它没有按 paper 的 Table 1 去跑完整 peak 列表。
4. 它外层只做了 5 次重复平均，不是 paper 里的 30 次 realization。
5. 它用了按图像切换的经验参数，而 paper 写的是统一参数 `d=13, D=11, mu=0.2, nu=1e-4, T=2`。
6. 它把其他方法都留在注释里，并没有真的提供 baseline 实现。

结论：**旧 `demo_onl.m` 只能当参考，不能当 paper 复现实验脚本。**

---

## 2.2 runner 必须改成与 paper 对齐

新 runner 只做下面几件事。

### A. 图像集合
目标图像应为：

- Spots
- Galaxy
- Cells
- Texture
- Fermi

实现策略：

1. 先检查 `images.mat` 里是否已经包含这 5 张图。
2. 如果 `Texture` 和 `Fermi` 已经在 `images.mat` 里，就直接复用。
3. 如果没有，就**不要重做 `images.mat`**，只额外新增一个简单 loader，把缺的图从 `data/` 目录单独读进来。

最小原则：**能不改 `images.mat` 就不改。**

### B. peak 列表
用于复现 Table 1 时，peak 直接按：

```matlab
peaks = [0.5 1 2 3 4 5];
```

注意：正文里对 peak 的文字描述和表格不完全一致。真正做复现时，以 **Table 1 / Table 4 里的 peak 列表**为准。

### C. 每张图先缩放到目标 peak
最小实现：

```matlab
u = peak * img / max(img(:));
```

说明：paper 说“scaled to the peak intensity level”，但没有在正文里把缩放代码写死。上面这个写法是最自然、也是最小的实现假设。

要求：

- 全程保持 `double`
- 不转 `uint8`
- 缩放后再 `poissrnd`

### D. 每个 image/peak 跑 30 次随机 realization

```matlab
num_trials = 30;
rng(0, 'twister');
```

每次都：

1. `z = poissrnd(u);`
2. `u1 = nlm_poisson004(z, 5, 6, 0.2);`
3. `u2 = imfilter(u1, fspecial('gaussian',5,1), 'symmetric');`
4. `nmise = fnmise(u, u2);`

最后统计：

```matlab
mean_nmise = mean(vals);
std_nmise  = std(vals);
```

### E. 第二步 Gaussian smoothing 固定为 paper 的 Step 2
不要再按图像切换 `nlwid/nlsi`。

统一写成：

```matlab
h = fspecial('gaussian', 5, 1);
out = imfilter(out1, h, 'symmetric');
```

这正对应 paper 的：

- `T = 2`
- 权重 `exp(-||x-x0||^2 / 2)`

### F. 去掉可视化副作用
把 `imshow(z,[])` 删掉或注释掉，保证脚本可以批处理、无界面运行。

---

## 2.3 对 `nlm_poisson004.c` 做最小论文对齐

这里是**唯一建议直接改源代码**的地方。

### 保留不动的部分

不要改这些：

- 镜像 padding
- `w` 的构造方式
- `avsqrt_image()` 这个预计算框架
- MEX 接口参数形式
- 主循环整体结构

### 只改三件事

#### 改动 1：显式加入 paper 里的 `nu = 1e-4`
在 `owf_denoising()` 内部增加常量：

```c
const double nu = 1e-4;
```

#### 改动 2：把当前 similarity 的实现改成更接近 paper Eq. (30)
当前代码是：

```c
dist2 = MAX((sqrt(dist2) - 1.414*v), 0);
```

这不是 paper 里写的
\((\text{weighted patch diff} - 2u_D(x_0))_+\)。

最小改法：

```c
double uD = v * v;   /* 因为 avimage[k] 存的是 sqrt(local mean) */
dist2 = MAX(dist2 - 2.0 * uD, 0.0);
```

也就是说：

- `dist2` 在进入这行前，仍表示“加权 patch 平方差”
- `uD` 用当前代码里已经有的局部均值估计
- 不要再对 `dist2` 先开方

#### 改动 3：把权重分母改成 paper 的 `H^2(x0)=mu*sqrt(uD)+nu`
当前代码是：

```c
wrho = exp(-dist2/(nlh*v));
```

最小改法：

```c
wrho = exp(-dist2 / (nlh * v + nu));
```

其中：

- `nlh` 就继续扮演 paper 里的 `mu`
- `v = sqrt(uD)` 已经由 `avsqrt_image()` 给出来了
- `nu` 直接 hardcode 为 `1e-4`

### 推荐的最小 patch 形状
在像素主循环里改成下面这个思路：

```c
v = avimage[k];
if (v == 0)
    denoisy[k] = 0;
else
{
    double uD = v * v;
    sum = 0;
    ar  = 0;

    for (i = 0, xp = x - nw; xp <= x + nw; xp++)
    {
        for (yp = y - nw; yp <= y + nw; yp++, i++)
        {
            adrp = yp * nxsy + xp;
            for (j = npnp, dist2 = 0., dd = dadr; j--; dd++)
            {
                dist  = synoisy[adr + *dd] - synoisy[adrp + *dd];
                dist2 += w[j] * dist * dist;
            }

            dist2 = MAX(dist2 - 2.0 * uD, 0.0);
            wrho  = exp(-dist2 / (nlh * v + nu));
            sum  += wrho;
            ar   += wrho * synoisy[adrp];
        }
    }
    denoisy[k] = MIN(MAX(ar / sum, 0), 255);
}
```

### 不要做的事

- 不要改函数签名
- 不要新加第五个输入参数去传 `nu`
- 不要改 memory layout
- 不要把 C 代码整个重写成 MATLAB 版
- 不要顺手重命名一堆变量

---

## 2.4 新增 `build_mex.m`

新增一个最小构建脚本：

```matlab
function build_mex()
    mex -O nlm_poisson004.c
end
```

如果编译器环境需要，可在本地手动先跑一次 `mex -setup`。

要求：**运行时只依赖当前机器重新编译出来的 MEX，不依赖 repo 自带的 `.mexw64/.dll`。**

---

## 3. 新 runner 的推荐结构

建议新增 `reproduce_table1_ours.m`，内容保持非常小。

## 3.1 伪代码

```matlab
function results = reproduce_table1_ours()
    build_mex();

    data = load_paper_images();   % 自己写一个很薄的 loader
    peaks = [0.5 1 2 3 4 5];
    num_trials = 30;

    D_half = 5;   % D = 11
    d_half = 6;   % d = 13
    mu = 0.2;

    h = fspecial('gaussian', 5, 1);
    rng(0, 'twister');

    rows = {};

    for p = 1:numel(peaks)
        peak = peaks(p);
        for k = 1:numel(data)
            img = data{k}.image;
            name = data{k}.name;
            u = peak * img / max(img(:));

            vals = zeros(num_trials, 1);
            for t = 1:num_trials
                z = poissrnd(u);
                u1 = nlm_poisson004(z, D_half, d_half, mu);
                u2 = imfilter(u1, h, 'symmetric');
                vals(t) = fnmise(u, u2);
            end

            rows(end+1,:) = {peak, name, mean(vals), std(vals)};
        end
    end

    results = cell2table(rows, ...);
    if ~exist('results','dir'), mkdir('results'); end
    writetable(results, fullfile('results', 'table1_ours.csv'));
    save(fullfile('results', 'table1_ours.mat'), 'results');
end
```

## 3.2 `load_paper_images.m` 的最小职责

它只做两件事：

1. 尝试从 `images.mat` 中取 `Spots/Galaxy/Cells/Texture/Fermi`
2. 如果缺 `Texture` 或 `Fermi`，就从 `data/` 下额外文件读取

不要把这个 loader 写成很复杂的数据下载器。

---

## 4. baseline 的处理原则

## 4.1 第一阶段完全不补 baseline

理由很简单：这个 repo 并没有把 paper 对比用到的下列方法完整打包进来：

- P-LET
- MV+B3
- MV-7/9
- P-NLM
- OWPF
- E+BM3D
- NLPCA
- FoEbin
- P4IP

`demo_onl.m` 里只有部分旧调用的注释痕迹，并不能直接跑出 Table 1 的完整对比。

所以最小改动方案里：

- **只复现 Ours 自己这一列**
- 然后再和 paper 的 Table 1 做人工对照

这是最现实、也最干净的最小方案。

## 4.2 如果后面一定要补 baseline

那是第二阶段，不应和这次最小改动混在一起。到时候另开 `baselines/` 目录，分别接入外部实现。

---

## 5. 结果校验

## 5.1 先做功能性校验

### MEX 校验

- `build_mex` 成功
- `nlm_poisson004(rand(32),5,6,0.2)` 返回与输入同尺寸矩阵

### 指标校验

- `fnmise(u,u) == 0`
- 脚本能完整跑完 `6 peaks x 5 images x 30 trials`

## 5.2 再做 paper 数值校验

不要要求逐位完全一致，但至少应满足：

1. 同一 peak 下，数值量级与 Table 1 接近。
2. 各图像上的相对强弱趋势基本一致。
3. 如果很多单元格偏差都很大，优先排查：
   - 图像集是否对了
   - peak 缩放是否对了
   - `dist2` 是否还保留了旧的 `sqrt(...) - 1.414*v`
   - 第二步 Gaussian 是否固定成了 `5x5, sigma=1`
   - trial 数是否真的等于 30

---

## 6. 建议输出文件

复现完成后，repo 里应当至少新增：

- `build_mex.m`
- `reproduce_table1_ours.m`
- `load_paper_images.m`（如果需要）
- `results/table1_ours.csv`
- `results/table1_ours.mat`

可选新增：

- `results/table4_ours.csv`（如果后续也算 PSNR/SSIM）

---

## 7. 不要改的东西

下面这些在最小改动复现里都不要碰：

- 不要把 repo 改成 Python 项目
- 不要把 `nlm_poisson004.c` 拆成很多文件
- 不要顺手重构变量名
- 不要顺手把 paper 里的所有 theorem 也代码化
- 不要把 baseline 一起塞进这次 patch
- 不要用现成 `.mexw64/.dll` 直接当最终结果

---

## 8. 如果没有 MATLAB 工具箱，才做的兜底改动

只有在缺工具箱时再补这些最小替代件：

- 没有 `poissrnd`：补一个很薄的 `poissrnd_local.m`
- 没有 `fspecial`：手写一个 5x5 Gaussian kernel
- 没有 `imfilter`：改成 `conv2(..., 'same')` 并自己做对称 padding

如果 MATLAB 环境正常，就不要为了“纯净”而额外重写这些。

---

## 9. 建议给 Codex 的任务说明（可直接复制）

```text
请在这个 repo 上做“最小改动复现”，要求如下：

1. 不换语言，不重构，不优化，只修复与论文不一致的地方。
2. 保留 fnmise.m 不变。
3. 保留 nlm_poisson004.c 的整体结构，只做最小论文对齐：
   - 在 owf_denoising() 里加入 const double nu = 1e-4;
   - 当前相似度实现是 MAX((sqrt(dist2)-1.414*v),0)，改成：
       uD = v*v;
       dist2 = MAX(dist2 - 2.0*uD, 0.0);
   - 当前权重是 exp(-dist2/(nlh*v))，改成：
       exp(-dist2 / (nlh*v + nu))
   - 不改函数签名，不加新输入参数。
4. 新增 build_mex.m，内容仅为 mex -O nlm_poisson004.c。
5. 不继续使用旧 demo_onl.m 作为论文复现脚本；新增 reproduce_table1_ours.m。
6. 新 runner 必须：
   - 使用图像：Spots, Galaxy, Cells, Texture, Fermi
   - 使用 peaks = [0.5 1 2 3 4 5]
   - 每个 image/peak 跑 30 次
   - 缩放方式：u = peak * img / max(img(:))
   - 调用：u1 = nlm_poisson004(z,5,6,0.2)
   - 第二步固定：fspecial('gaussian',5,1) + imfilter(...,'symmetric')
   - 指标：fnmise(u, u2)
   - 输出 results/table1_ours.csv 和 results/table1_ours.mat
   - 注释或删除 imshow
7. 如果 images.mat 里缺 Texture/Fermi，则不要重做 images.mat，只新增一个很薄的 load_paper_images.m，从 data/ 读取缺失图片。
8. 不要在这次 patch 里加入 baseline 实现。
9. 最终给出简短运行说明。
```

---

## 10. 人工校验用：Table 1 中 Ours 的目标值（手工抄录）

下表只用于复现后快速比对趋势，不要求逐位完全一样。

| Peak | Spots | Galaxy | Cells | Texture | Fermi | Average |
|---|---|---|---|---|---|---|
| 0.5 | 0.0149±0.0014 | 0.0075±0.0004 | 0.0105±0.0004 | 0.0095±0.0003 | 0.0087±0.0033 | 0.0102±0.0012 |
| 1 | 0.0243±0.0010 | 0.0128±0.0004 | 0.0174±0.0005 | 0.0155±0.0003 | 0.0133±0.0019 | 0.0157±0.0008 |
| 2 | 0.0257±0.0012 | 0.0183±0.0004 | 0.0306±0.0006 | 0.0260±0.0004 | 0.0315±0.0012 | 0.0264±0.0008 |
| 3 | 0.0311±0.0009 | 0.0204±0.0004 | 0.0324±0.0006 | 0.0356±0.0005 | 0.0335±0.0013 | 0.0306±0.0007 |
| 4 | 0.0322±0.0012 | 0.0262±0.0005 | 0.0449±0.0008 | 0.0448±0.0006 | 0.0452±0.0015 | 0.0387±0.0009 |
| 5 | 0.0338±0.0011 | 0.0282±0.0005 | 0.0411±0.0007 | 0.0531±0.0006 | 0.0455±0.0024 | 0.0403±0.0011 |

如果你的结果整体趋势差得很远，先不要调参，先检查：

- 图像是否换成了 paper 的 5 张图
- peak 是否用了 `[0.5 1 2 3 4 5]`
- trial 是否真的是 30
- C 代码里是否还保留旧的 `sqrt(dist2)-1.414*v`
- 第二步 Gaussian 是否固定成 `5x5, sigma=1`

