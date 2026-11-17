# Paired FASTQ Re-pairing Service

测序湖入湖工具：接收顺序被打乱的双端未压缩 FASTQ（R1/R2），按读段号重新
配对，按 R1 顺序交付两份文件，原文逐字节保留。FastAPI + SQLite 实现，单文件
上限 256 MiB，进程内存上限 256 MiB，可处理总量超过内存的输入。

服务还提供**光学重复核查**：对成功任务的结果文件逐对核对，按像素距离与双端
序列聚类，给出待去重读号清单及其保留读号（见下文「光学重复核查」）。

## 配对规则

- 读段号取头行首个空白前的字段，去掉开头 `@`，再去掉末尾的 `/1` 或 `/2`
  （例如 `@read-7/1 sample=x` → `read-7`）。
- 输出按 **R1 的原始顺序**排列；R2 按相同读段号对齐到同一位置。
- 四行记录原文（含头行描述、换行风格）原样写入结果，不做任何改写。

## 失败即整次失败，不开放部分结果

任一情况发生，任务整体标记为 `failed`/`cancelled`，已写入的暂存全部删除，
下载接口返回 409：

- 同一文件内读段号重复（含 `/1`、`/2` 归一化后重复）
- 孤端（某读段号只出现在一侧）
- 格式错误：非四行式、头行不以 `@` 开头、第三行不以 `+` 开头、空文件
- 序列与质量字符串长度不等

## 内存与大输入

- 上传经 `SpooledTemporaryFile` 落盘（超过 1 MiB 即写临时文件），分块复制，
  单文件流式校验 256 MiB 上限，超限返回 413。
- 解析按 512 条/批写入 SQLite（WAL、小 page cache、临时表走文件
  `temp_store=FILE`、关闭 mmap），任意时刻内存里只有一个批次。
- 配对发布使用 `JOIN ... ORDER BY a.ord` 流式游标，分批 `fetchmany` 写出。
- 因此输入总量可以远超 256 MiB；自测以约 393 MiB 总量（两份各 196.5 MiB、
  各 200 万条）验证，峰值 RSS 约 19 MiB（见下文自测结果）。

## 任务与取消语义

- `POST /tasks` 上传后立即返回 `202`，后台守护线程处理。
- `GET /tasks/{id}` 查询状态、阶段（upload/parse/validate/publish/done）与
  0–1 进度。
- `POST /tasks/{id}/cancel`：置取消标记后**立即返回**，不等待配对结束；
  工作线程在下一个批次边界中止。取消与完成之间用条件 UPDATE 仲裁：
  - 取消先到 → 任务最终 `cancelled`，结果永不发布，暂存目录删除；
  - 完成先到 → 取消请求得到 `409`，明确告知任务已完成、不可取消，结果保留。
- 取消后可立刻创建新任务；任务目录以 UUID 隔离，清理只删除本任务目录，
  不会误删新任务文件。
- 服务重启（启动钩子）：所有非终态任务标记 `failed`（原因
  `interrupted by service restart`）并清理暂存与记录；`succeeded` 任务的
  结果文件保留，仍可下载。

## API

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/tasks` | multipart 字段 `r1`、`r2`，返回任务对象（202） |
| GET  | `/tasks/{id}` | 状态/进度/消息/下载链接 |
| POST | `/tasks/{id}/cancel` | 请求取消（立即返回） |
| GET  | `/tasks/{id}/download/r1` | 下载配对后的 R1 |
| GET  | `/tasks/{id}/download/r2` | 下载配对后的 R2 |
| GET  | `/health` | 健康检查 |
| POST | `/optical-checks` | 对成功任务做光学重复核查（同步返回完整清单） |
| GET  | `/optical-checks/{id}` | 取回已生成的核查报告 |

## 光学重复核查

输入为某次**成功**配对任务的 `task_id` 与像素距离阈值（1–10 的整数）。核查
直接读取该任务已发布的结果文件，按原配对规则逐行核对 R1/R2：

- 两条读对若 **同仪器、同批次（run）、同流动槽、同泳道、同瓦片** 且
  **双端序列完全相同**，则视为同一候选组；
- 组内两条读对的坐标 `(X, Y)` 欧氏距离 `≤ 阈值` 时连边；连通的读对归为一簇
  （传递闭包，链式相连也算同簇）；
- 每簇保留 **双端 Phred+33 质量总分最高** 的读对；同分保留在 R1 中靠前者；
- 报告按 **R1 顺序** 给出全部待去重读号及其保留读号，以及总对数与重复对数。
- 读段号须为 `INST:RUN:FLOWCELL:LANE:TILE:X:Y`（归一化规则与配对一致：去掉
  开头 `@`、去掉末尾 `/1`、`/2`；LANE/TILE/X/Y 为十进制非负整数，LANE/TILE
  至少为 1）。

```bash
curl -X POST http://localhost:8000/optical-checks \
  -H 'Content-Type: application/json' \
  -d '{"task_id": "<成功任务号>", "distance_threshold": 2}'
```

响应示例（`201`）：

```json
{
  "id": "…",
  "task_id": "…",
  "state": "ready",
  "distance_threshold": 2,
  "total_pairs": 4,
  "duplicate_pairs": 2,
  "duplicates": [
    {"read_id": "A:1:FC:1:100:10:10", "keep_read_id": "A:1:FC:1:100:10:11"},
    {"read_id": "A:1:FC:1:100:13:14", "keep_read_id": "A:1:FC:1:100:10:11"}
  ],
  "error": "",
  "created_at": "…",
  "updated_at": "…"
}
```

失败语义与配对一致——**整次失败，不返回部分清单，也不改动原件**：

- 任务不存在 → `404`；任务非 `succeeded` 或结果文件缺失 → `409`；
- 阈值不在 1–10、`task_id` 类型错误 → `422`（请求校验，未启动核查）；
- 记录损坏（非四行、头行/第三行前缀错误、截断、R1/R2 条数或读号不一致、
  重复读号）、坐标非法、读段号不符合上述格式、质量字符不在 ASCII 33–126、
  单端超过 500 碱基、总对数超过 50,000 → `422`，`detail` 说明原因并给出
  `check_id`；数据库里该次核查记为 `failed`，`GET /optical-checks/{id}` 可查；
- 同时只允许一个核查在跑（内存上限约束）；第二个请求得到 `429`。

结果文件只读，核查前后内容与修改时间均不变（自测用 SHA-256 前后比对）。

内存：核查只在内存里保留每对的轻量元数据（读号、坐标、质量总分，序列只留
16 字节摘要），不保留序列与质量原文。自测最坏情况——50,000 对、每端 500
碱基、50,000 个不同序列组、阈值 10——单独调用分析器峰值 RSS 约 50 MiB；
经完整 uvicorn 服务实测进程峰值 RSS 约 104 MiB，均低于 256 MiB；耗时约
23–26 秒。

## 运行（Docker Compose，内存上限 256 MiB）

```bash
docker compose up --build
# 打开 http://localhost:8000/docs
```

`docker-compose.yml` 中 `mem_limit: 256m` 即交付要求的内存上限；数据保存在
命名卷 `fastq-data`（容器内 `/data`，含 SQLite 与各任务目录），重启后成功
结果仍可下载。

## 本地直接运行（开发/自测）

```bash
python3 -m venv --without-pip .venv
.venv/bin/python tools/get-pip.py
.venv/bin/pip install -r requirements.txt
FASTQ_DATA_DIR=./data .venv/bin/uvicorn app.main:app --port 8000
```

## 演示脚本

```bash
./scripts/demo.sh                 # 需要一个运行在 localhost:8000 的服务
BASE_URL=http://localhost:8000 ./scripts/demo.sh
```

脚本会生成一对乱序 FASTQ，演示上传 → 轮询进度 → 下载校验，并额外演示一次
取消、一次错误输入、一次光学重复核查（含一个光学重复簇）以及一次因读段号
不符合 `INST:RUN:FLOWCELL:LANE:TILE:X:Y` 而整次失败的核查。

## 自测

```bash
# 短测试（配对/取消/重启 + 光学重复核查 + 取消与完成仲裁），共 102 项检查
bash scripts/run_tests.sh
# 或分别运行
.venv/bin/python -I tests/test_api.py          # 30 项：配对、整次失败、取消、重启恢复
.venv/bin/python -I tests/test_optical.py      # 63 项：光学重复核查（聚类/保留规则/整次失败/边界）
.venv/bin/python -I tests/test_cancel_race.py  # 9 项：取消/完成条件 UPDATE 仲裁

# 内存压力测试：两份各 196.5 MiB（合计约 393 MiB，超过 256 MiB 内存上限），
# 子进程实测峰值 RSS 约 19 MiB，并逐条校验配对顺序与原文
.venv/bin/python -I tests/test_memory.py 2000000
```

## 目录结构

```
app/
  config.py    # 上限、路径、批大小等配置（环境变量可覆盖）
  db.py        # SQLite schema 与任务/记录/核查报告访问
  fastq.py     # 流式解析、校验、JOIN 配对发布
  optical.py   # 光学重复核查：结果文件核对、聚类、保留规则
  worker.py    # 后台线程、取消仲裁、启动恢复、任务目录清理
  main.py      # FastAPI 路由与有界流式上传、光学核查接口
tests/         # 功能自测、光学核查自测与内存压力测试
scripts/demo.sh
Dockerfile, docker-compose.yml, requirements.txt
```

## 遗留项 / 边界说明

- 仅接受未压缩 FASTQ；`.gz` 等压缩格式不在范围内。
- 服务为单实例 + 单 SQLite 文件设计；多实例横向扩展需要换用共享存储/队列。
- 没有内置任务 TTL 清理，成功结果会一直保留在数据卷中，需外部运维定期清理。
- 光学核查同步执行：最坏 50,000 对约 23 秒，期间该请求占用一个工作线程；
  同一时刻只跑一个核查（第二个请求 `429`）。若需要后台化/进度条可改为任务式。
- 光学核查要求结果文件成对对齐（本服务发布的结果天然满足）；若结果被外部
  改动导致顺序错位，按整次失败处理。
- 开发沙箱内未安装 Docker，`docker compose up` 未实跑；compose/Dockerfile 仅做
  了 YAML 与结构校验，应用本身已用等价的本地 uvicorn 命令完整验证（含
  `/docs`、上传、取消、下载、重启恢复、光学核查）。在装有 Docker 的环境执行
  `docker compose up --build` 即可，`mem_limit: 256m` 即内存上限。
