# 只上传这个新文件夹

先确认已停止运行本目录的仿真，再在 PowerShell 中执行：

```powershell
Set-Location 'G:\2026Huawei\Code\Math_Modeling_2026_B'
git status --short
git add -- 'Code_Q2_PhaseFix/Q2_Part2_Step_Loading_Optimization_Verified'
git diff --cached --name-only
```

确认暂存清单**只有这个新文件夹**，没有原模型目录、build_grid.m 或队友文件；否则先停止，不要提交。确认无误后：

```powershell
git commit -m "Add audited step loading optimization package"
git -c http.proxy=http://127.0.0.1:7890 -c http.sslBackend=openssl push origin main
```

代理仅在本机对应服务启用时使用；否则使用 `git push origin main`。若推送提示远端更新，不要使用force；先检查本地改动，再安全同步远端。

整个文件夹可独立使用。`.gitignore`只忽略可续跑缓存、编辑器备份和日志；真实验证MAT、参考曲线MAT和报告可以随包提交。`MODEL_SHA256.txt`记录本地交付时模型文件字节校验值，Git行尾转换可能改变字节哈希，不代表方程变化。

Claude原文件夹和既有修改未被本修复覆盖；本任务没有自动运行上述提交/推送命令。
