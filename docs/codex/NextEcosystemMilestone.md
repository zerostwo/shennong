# Shennong 下一阶段执行任务书

创建日期：2026-07-28

状态：可在新的 Codex 窗口中执行

## 目标

下一阶段只做一件主线任务：

> 使用已发布的不可变镜像部署 ShennongDB、Shennong Runtime 和 Shennong
> OS，并用真实 PBMC3K 数据完成第一次五仓端到端分析、结果提升、回读和
> lineage 验证。

本阶段完成前，不得宣称 scRNA-seq 已实现生态级端到端支持，也不得把通用
文件传输能力描述为空间、CITE-seq、scATAC-seq 或通用 bulk omics 支持。

## 必读入口

新窗口开始后依次阅读：

1. `/home/duansq/dev/packages/shennong/AGENTS.md`
2. `/home/duansq/dev/packages/shennong/docs/codex/Ecosystem.md`
3. `/home/duansq/dev/packages/shennong/docs/codex/ecosystem-lock.json`
4. 五个仓库各自的 `AGENTS.md`
5. 现有部署、Compose、迁移和恢复文档

`/home/duansq/dev/services/shennong.one` 不是 Git 仓库。所有 Git
操作、验证、提交和推送必须在五个仓库中分别完成。

## 当前锁定基线

开始执行前必须重新核对，不能仅凭本文件假定它们仍然有效。

| 组件 | 锁定 revision |
|---|---|
| Shennong 实现 | `e58f7fbfea1ea02cebc07efddcabb66d6e8bb48a` |
| Shennong 协调仓库 HEAD | `51c7ae490fa0ebba6a59369b944871a1e8c1204c` |
| Shennong Runtime 内的 Shennong pin | `c1d958db3319f635ff5d6f9ad484a208774a4a39` |
| ShennongData | `17f0f0e87dd8ad2a3751dd11c58c8aa43823aa69` |
| Shennong OS | `e3f421bcba70e85a82596634178f24f3270df621` |
| Shennong Runtime | `1649f1da43b25d79e8d820bd1ef093cf49a77114` |
| ShennongDB | `24445c3b08bb708bcdc2d1cf2eccba1735815633` |

部署必须使用以下 OCI index digest，不得使用浮动的 `latest`：

| 服务 | 不可变镜像 |
|---|---|
| OS | `zerostwo/shennong-os@sha256:718b5be9e42c4dda71c57dec78c2fdcf0533a6a40eb94edd5cee498f191c9615` |
| Runtime | `zerostwo/shennong-runtime@sha256:6506ce58451b15cd100ac23e14524f75aae2b4c120bc7852a4e6bd8ee1e0226b` |
| DB | `zerostwo/shennong-db@sha256:0c4220718d9fe5684dcba963c1b9abad9d753b4a4b86177067b5e5fef6b4e5f3` |

Runtime 工具链 SHA-256：

`ea92fd5c3c8cf2b4da5b5e409cde94179f6c358e0e97eac178be151fe5eda6bf`

当前兼容锁仍是 `candidate`，所有服务的 `deployed` 均为 `false`。

## 不可破坏的架构边界

- 用户和 R 客户端路径是
  `ShennongData -> Shennong OS -> ShennongDB`。
- ShennongData 不得获得 DB admin key。
- Runtime 不得获得通用 DB 凭据。
- OS 是用户、Project、授权、Job 和平台执行 provenance 的权威。
- DB 是不可变 Resource、revision、Artifact bytes 和 lineage 的权威。
- Runtime 是隔离执行、实际输出文件、大小和摘要的权威。
- producer Result Bundle 的执行字段只是可选一致性声明，不能创建平台
  Project、actor、Run、Job、Artifact 或 Activity 身份。
- projectless 调用仅允许公开目录发现；inspect、query、resolve、download
  必须提供 Project UUID。
- 任何凭据、PAT、cookie、admin key、预签名 URL 或用户数据都不得写入
  Git、日志、兼容锁或最终报告。

## 阶段 0：只读预检

在任何部署或迁移前：

1. 在五个仓库分别记录：
   - 当前分支与 HEAD；
   - `HEAD...@{upstream}` divergence；
   - tracked/untracked/ignored 变化；
   - 当前 CI、发布镜像和 live deployment 状态。
2. 保留 ShennongData 中既有的未跟踪文件
   `ShennongData_R_client_architecture_design.md`，除非用户明确要求处理。
3. 核对三个 OCI index、amd64 manifest 和
   `org.opencontainers.image.revision`。
4. 查明唯一的目标部署环境、现有 Compose/服务编排、数据库位置、备份与
   rollback 方法。
5. 如果目标环境不唯一、迁移前无法备份、锁定 revision 已漂移，或存在
   与本任务重叠的未提交修改，停止部署并向用户说明。
6. 记录当前服务版本、健康状态、端口和数据库 migration 状态，作为回滚
   前基线。

不要因为“CI 绿色”或“镜像已发布”跳过这一阶段。

## 阶段 1：按不可变 digest 部署

推荐顺序：

1. **ShennongDB**
   - 部署锁定 digest；
   - 备份后执行已有 migration；
   - 验证健康、API/OpenAPI、storage backend、当前 revision；
   - 验证 local/S3 模式下精确 Artifact metadata 和 content 读取。
2. **Shennong Runtime**
   - 部署锁定 digest；
   - 验证 `/v1/health` 和 `/v1/info`；
   - `/v1/info` 必须返回锁定的 Runtime revision、两个 R 包 commit、
     package version 和 toolchain SHA-256。
3. **Shennong OS**
   - 部署锁定 digest；
   - 执行 OS 自身 migration；
   - 验证 `/healthz`、Web、Agent、DB、Runtime 连接；
   - 确认 OS 的 Runtime admission 使用完整工具链锁，而非只检查 SemVer。

每一步都必须：

- 使用现有部署/恢复脚本，不临时发明第二套生产拓扑；
- 保存脱敏后的命令和响应证据；
- 验证 rollback 路径；
- 失败时停止后续组件升级，不能把部分部署标记为完成。

## 阶段 2：部署后安全与接口 smoke test

至少验证以下正向和反向路径：

### PAT 与 Project

1. 创建最小权限 PAT，确认明文只出现一次且响应为 `no-store`。
2. 使用 PAT 完成同 Project 的公开发现和受保护查询。
3. 无 Project 的 inspect/query/resolve/download 必须失败。
4. 未绑定或其他 Project 的 private Resource 必须 fail closed。
5. 撤销 PAT 后，同一请求必须立即失败。

### 公开目录

- 匿名目录只能返回白名单字段；
- 不得泄露 `spec`、provenance、owner、checksum、凭据或 storage/source
  locator；
- 匿名详情和 child 请求必须失败。

### Artifact

- DB 的 `X-Shennong-No-Redirect: 1` 只对经过认证、无用户身份的内部
  service administrator 生效；
- 用户、错误凭据和未认证调用不能获得内部 byte proxy；
- OS 必须重新校验大小和 SHA-256；
- Range 请求按当前合同返回 416，OS/BFF 不宣告部分读取能力；
- 零字节 Artifact 的 `Content-Length: 0` 和空内容 SHA-256 必须正确。

所有反向测试也必须作为发布证据保留。

## 阶段 3：PBMC3K 五仓端到端 fixture

目标路径：

```text
ShennongDB immutable PBMC3K Resource/Artifact
  -> Shennong OS Project/PAT authorization
  -> ShennongData discovery, download and Data Bundle v1
  -> Shennong Runtime exact toolchain admission
  -> Shennong scRNA-seq analysis and result validation
  -> Result Bundle v1
  -> OS authoritative promotion
  -> DB immutable result Artifact and Activity lineage
  -> Project-scoped discovery and exact result readback
```

### 3.1 输入权威与物化

1. 使用仓库已有的真实 PBMC3K/10x MEX provider，不创建伪造的
   “PBMC-like” CSV。
2. 解析并记录精确 Resource revision、Artifact UUID、字节大小和
   SHA-256。
3. 通过 OS-issued PAT 和 Project UUID 发现输入。
4. 使用已安装的 ShennongData 生成
   `shennong.dev/data-bundle/v1`。
5. 验证：
   - feature/cell 标识唯一且稳定；
   - count 非负，缺失 feature 与结构零没有混淆；
   - matrix、cell metadata、feature metadata 和 immutable provenance
     一致；
   - 下载内容的大小和 SHA-256 与 DB 完全一致。

### 3.2 Runtime 与 Shennong 执行

1. OS 创建绑定 Project、actor、Job、输入 Artifact 和工具链锁的执行。
2. Runtime admission 必须精确匹配锁定 schema、版本、commit 和
   toolchain SHA-256。
3. 运行一个确定性的最小 scRNA-seq 流程，至少覆盖：
   - Seurat 初始化；
   - QC 和 normalization；
   - dimensional reduction / clustering；
   - 一个可验证的结果表或图形 Artifact。
4. 固定随机种子并记录实际方法、参数、包版本和 Runtime revision。
5. 使用 `sn_validate_result()` 验证分析结果。
6. 使用 `sn_build_result_bundle()` /
   `sn_validate_result_bundle()` /
   `sn_export_result_bundle()` 生成真实的
   `shennong.dev/analysis-result-bundle/v1`，不能手写一个看似正确的
   JSON 替代 R producer fixture。

### 3.3 科学断言

断言应从 PBMC3K manifest 和实际对象推导，不要凭记忆硬编码维度：

- 输入 features/cells 与 manifest 一致；
- counts 非负且数据不是空矩阵；
- QC 指标有限且覆盖所有 cells；
- PCA/UMAP 或选定 embedding 的行与 cell IDs 对齐且数值有限；
- clustering 至少产生一个非空 cluster assignment，且没有丢失 cells；
- canonical analysis result schema 有效；
- 输出表、图或序列化对象可以独立重新读取；
- 同一固定输入、参数和 seed 的关键结果可重现。

如果只能证明软件执行成功而不能证明输入与结果科学有效，不得通过本阶段。

### 3.4 提升、lineage 与回读

1. Runtime 从实际文件计算 canonical output identifier、size 和 SHA-256。
2. OS 以 Job/Project/actor/manifest/DB state 为权威，验证 producer claims。
3. OS 上传结果 bytes，按 deterministic filename 精确回读并复核：
   - actor；
   - Project；
   - filename；
   - content type；
   - size；
   - SHA-256；
   - storage URI。
4. DB 创建不可变 result Resource/Artifact 及 Activity input/output
   lineage。
5. 只有在提升、摘要回读和 lineage 全部成功后，Job/plan step 才可完成。
6. 再通过 Project-scoped discovery 找到结果，并下载、验证、重新读取。

## 阶段 4：证据与收口

创建本次执行报告，例如：

`docs/codex/runs/20260728-pbmc3k-five-repo-e2e.md`

报告必须包括：

- 五仓执行时的 branch、HEAD、upstream 和 dirty state；
- 三个实际部署镜像 digest 与 live revision 响应；
- migration、备份和 rollback 证据；
- PAT 正向/反向测试；
- 输入 Resource/revision/Artifact 与 SHA-256；
- Data Bundle、Runtime toolchain 和 Result Bundle schema；
- Shennong 方法、参数、seed 和科学断言；
- 结果 Resource/Artifact/Activity lineage；
- 精确回读的 size/SHA-256；
- 所有测试、CI、commit、push 和部署结果；
- 脱敏后的失败与恢复记录。

随后按实际证据更新：

- `docs/codex/Ecosystem.md`
- `docs/codex/Status.md`
- `docs/codex/Decisions.md`
- `docs/codex/ecosystem-lock.json`

只有三个 live 服务都报告锁定 revision 时，才可以将对应
`deployed` 改为 `true`。只有完整 PBMC3K 路径和科学断言全部通过时，才可以
将 `scrna_seq.all_five_repositories_end_to_end` 改为 `true`。

每个仓库独立执行测试、更新 changelog、使用 Conventional Commit、推送并
记录远端 CI。不得从 `shennong.one` 聚合目录执行 Git 操作，不得 force
push，不得覆盖无关工作。

## 明确不在本阶段内

- 稳定 cursor pagination；
- bulk RNA 五仓 fixture；
- 空间 companion/image/FOV/segmentation 合同；
- CITE-seq RNA+ADT 五仓 fixture；
- scATAC 分析 API；
- 通用 proteomics、metabolomics 或 epigenomics 平台支持。

这些进入 PBMC3K 闭环之后的后续里程碑。发现阻断 PBMC3K 的真实缺口时，可以
做最小修复，但不得顺便扩展成新模态工程。

## Definition of Done

只有同时满足以下条件才算完成：

- [ ] 五仓预检完成且无未解释的 revision/dirty-state 漂移；
- [ ] DB、Runtime、OS 按锁定 digest 部署；
- [ ] live health/info/revision 与兼容锁一致；
- [ ] migration、备份和 rollback 证据完整；
- [ ] PAT、Project isolation、公开白名单和 no-redirect 反向测试通过；
- [ ] 真实 PBMC3K 经 OS 和 ShennongData 生成 Data Bundle v1；
- [ ] Runtime 使用精确工具链执行真实 Shennong scRNA 流程；
- [ ] 真实 R producer Result Bundle v1 通过验证；
- [ ] 结果 bytes、size、SHA-256、Resource/Artifact/Activity lineage 完整；
- [ ] 结果能够通过 Project-scoped 路径发现、下载和重新读取；
- [ ] 科学断言通过；
- [ ] 文档、兼容锁、changelog、commit、push 和远端 CI 完成；
- [ ] 没有凭据或隐私数据进入 Git 或日志。

## 新窗口启动指令

在新的 Codex 窗口中发送：

> 请先完整阅读
> `/home/duansq/dev/packages/shennong/AGENTS.md` 和
> `/home/duansq/dev/packages/shennong/docs/codex/NextEcosystemMilestone.md`，
> 然后严格按任务书执行。先做五仓只读预检和部署目标确认；保留所有无关
> 修改。完成不可变 digest 部署、PBMC3K 五仓端到端、科学断言、lineage
> 回读、文档、提交、推送和远端 CI 验收。在缺少唯一部署目标、备份或必要
> 权限时，停在首次破坏性操作前向我报告，不要猜测。
