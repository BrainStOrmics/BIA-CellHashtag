# Backends

[Core capabilities](/oss/python/deepagents/harness)

# Backends
Copy page

Choose and configure filesystem backends for Deep Agents. You can specify routes to different backends, implement virtual filesystems, and enforce policies.Copy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Deep Agents expose a filesystem surface to the agent via tools like `ls`, `read_file`, `write_file`, `edit_file`, `glob`, and `grep`. These tools operate through a pluggable backend. The `read_file` tool natively supports image files (`.png`, `.jpg`, `.jpeg`, `.gif`, `.webp`) across all backends, returning them as multimodal content blocks.
Sandboxes and the [`LocalShellBackend`](https://reference.langchain.com/python/deepagents/backends/local_shell/LocalShellBackend) also provide an `execute` tool.
This page explains how to:

{' ' * (self.list_depth - 1)}- 
[choose a backend](#specify-a-backend),

{' ' * (self.list_depth - 1)}- 
[route different paths to different backends](#route-to-different-backends),

{' ' * (self.list_depth - 1)}- 
[implement your own virtual filesystem](#use-a-virtual-filesystem) (e.g., S3 or Postgres),

{' ' * (self.list_depth - 1)}- 
[set permissions](#permissions) on filesystem access,

{' ' * (self.list_depth - 1)}- 
[comply with the backend protocol](#protocol-reference),

## [​](#quickstart)Quickstart

Here are a few prebuilt filesystem backends that you can quickly use with your deep agent:

| 
 | Built-in backend | Description
 | [Default](#statebackend) | `agent = create_deep_agent(model="google_genai:gemini-3.5-flash")` 
 Thread-scoped. The default filesystem backend for an agent is stored in `langgraph` state. Files persist across turns within a thread (via your checkpointer) and are not shared across threads.
 | [Local filesystem persistence](#filesystembackend-local-disk) | `agent = create_deep_agent(model="google_genai:gemini-3.5-flash", backend=FilesystemBackend(root_dir="/Users/nh/Desktop/"))` 
This gives the deep agent access to your local machine’s filesystem. You can specify the root directory that the agent has access to. Note that any provided `root_dir` must be an absolute path. Typically, wrap in a [CompositeBackend](#compositebackend-router) to keep internal agent data (offloaded tool results, conversation history) separate from your project files.
 | [Durable store (LangGraph store)](#storebackend-langgraph-store) | `agent = create_deep_agent(model="google_genai:gemini-3.5-flash", backend=StoreBackend())` 
This gives the agent access to long-term storage that is *persisted across threads*. This is great for storing longer term memories or instructions that are applicable to the agent over multiple executions.
 | [Context Hub](#contexthubbackend) | `agent = create_deep_agent(model="google_genai:gemini-3.5-flash", backend=ContextHubBackend("my-agent"))` 
Stores files durably in a LangSmith Hub repo, without provisioning a separate LangGraph store.
 | [Sandbox](/oss/python/deepagents/sandboxes) | `agent = create_deep_agent(model="google_genai:gemini-3.5-flash", backend=sandbox)` 
Execute code in isolated environments. Sandboxes provide filesystem tools plus the `execute` tool for running shell commands. Choose from Modal, Daytona, Deno, or local VFS.
 | [Local shell](#localshellbackend-local-shell) | `agent = create_deep_agent(model="google_genai:gemini-3.5-flash", backend=LocalShellBackend(root_dir=".", env={"PATH": "/usr/bin:/bin"}))` 
Filesystem and shell execution directly on the host. No isolation—use only in controlled development environments. See [security considerations](#localshellbackend-local-shell) below.
 | [Composite](#compositebackend-router) | Thread-scoped by default, `/memories/` persisted across threads. The Composite backend is maximally flexible. You can specify different routes in the filesystem to point towards different backends. See Composite routing below for a ready-to-paste example.

## [​](#built-in-backends)Built-in backends

### [​](#statebackend)StateBackend

```
from deepagents import create_deep_agent
from deepagents.backends import StateBackend

# By default we provide a StateBackend
agent = create_deep_agent(model="google_genai:gemini-3.5-flash")

# Under the hood, it looks like
agent2 = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=StateBackend(),
)

```

**How it works:**

{' ' * (self.list_depth - 1)}- Stores files in LangGraph agent state for the current thread via [`StateBackend`](https://reference.langchain.com/python/deepagents/backends/state/StateBackend).

{' ' * (self.list_depth - 1)}- Persists across multiple agent turns on the same thread via checkpoints. Files are not shared across threads.

Designed to be used from within a graph. Calling backend methods (e.g., `state_backend.upload_files(...)`) outside of a graph run won’t take effect until the graph executes.
**Best for:**

{' ' * (self.list_depth - 1)}- A scratch pad for the agent to write intermediate results.

{' ' * (self.list_depth - 1)}- Automatic eviction of large tool outputs which the agent can then read back in piece by piece.

Note that this backend is shared between the supervisor agent and subagents, and any files a subagent writes will remain in the LangGraph agent state
even after that subagent’s execution is complete. Those files will continue to be available to the supervisor agent and other subagents.

### [​](#filesystembackend-local-disk)FilesystemBackend (local disk)

[`FilesystemBackend`](https://reference.langchain.com/python/deepagents/backends/filesystem/FilesystemBackend) reads and writes real files under a configurable root directory.
This backend grants agents direct filesystem read/write access.
Use with caution and only in appropriate environments.**Appropriate use cases:**

{' ' * (self.list_depth - 1)}- Local development CLIs (coding assistants, development tools)

{' ' * (self.list_depth - 1)}- CI/CD pipelines (see security considerations below)
**Inappropriate use cases:**

{' ' * (self.list_depth - 1)}- Web servers or HTTP APIs - use `StateBackend`, `StoreBackend`, or a [sandbox backend](/oss/python/deepagents/sandboxes) instead
**Security risks:**

{' ' * (self.list_depth - 1)}- Agents can read any accessible file, including secrets (API keys, credentials, `.env` files)

{' ' * (self.list_depth - 1)}- Combined with network tools, secrets may be exfiltrated via SSRF attacks

{' ' * (self.list_depth - 1)}- File modifications are permanent and irreversible
**Recommended safeguards:**

{' ' * (self.list_depth - 1)}- Enable [Human-in-the-Loop (HITL) middleware](/oss/python/deepagents/human-in-the-loop) to review sensitive operations.

{' ' * (self.list_depth - 1)}- Exclude secrets from accessible filesystem paths (especially in CI/CD).

{' ' * (self.list_depth - 1)}- Use a [sandbox backend](/oss/python/deepagents/sandboxes) for production environments requiring filesystem interaction.

{' ' * (self.list_depth - 1)}- **Always** use `virtual_mode=True` with `root_dir` to enable path-based access restrictions (blocks `..`, `~`, and absolute paths outside root).
Note that the default (`virtual_mode=False`) provides no security even with `root_dir` set.

```
from deepagents import create_deep_agent
from deepagents.backends import FilesystemBackend

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=FilesystemBackend(root_dir=".", virtual_mode=True),
)

```

**How it works:**

{' ' * (self.list_depth - 1)}- Reads/writes real files under a configurable `root_dir`.

{' ' * (self.list_depth - 1)}- You can optionally set `virtual_mode=True` to sandbox and normalize paths under `root_dir`.

{' ' * (self.list_depth - 1)}- Uses secure path resolution, prevents unsafe symlink traversal when possible, can use ripgrep for fast `grep`.

**Best for:**

{' ' * (self.list_depth - 1)}- Local projects on your machine

{' ' * (self.list_depth - 1)}- CI sandboxes

{' ' * (self.list_depth - 1)}- Mounted persistent volumes

**Wrap `FilesystemBackend` in a `CompositeBackend`** for most use cases. Deep Agents automatically write internal data to the backend, including offloaded large tool results (under `/large_tool_results/`) and conversation history (under `/conversation_history/`). When you use `FilesystemBackend` alone, these internal files are written to real disk under `root_dir`, mixing agent artifacts with your project files.Use a `CompositeBackend` to route your project directory to `FilesystemBackend` while keeping internal paths in ephemeral `StateBackend` storage:
```
from deepagents import create_deep_agent
from deepagents.backends import CompositeBackend, StateBackend, FilesystemBackend

agent = create_deep_agent(
 backend=CompositeBackend(
 default=StateBackend(),
 routes={
 "/workspace/": FilesystemBackend(root_dir="/path/to/project", virtual_mode=True),
 },
 )
)

```
This way, agent reads and writes under `/workspace/` go to real disk, while offloaded tool results and other internal data stay in ephemeral state. See [Route to different backends](#route-to-different-backends) for more routing patterns.

### [​](#localshellbackend-local-shell)LocalShellBackend (local shell)

This backend grants agents direct filesystem read/write access **and** unrestricted shell execution on your host.
Use with extreme caution and only in appropriate environments.**Appropriate use cases:**

{' ' * (self.list_depth - 1)}- Local development CLIs (coding assistants, development tools)

{' ' * (self.list_depth - 1)}- Personal development environments where you trust the agent’s code

{' ' * (self.list_depth - 1)}- CI/CD pipelines with proper secret management
**Inappropriate use cases:**

{' ' * (self.list_depth - 1)}- Production environments (such as web servers, APIs, multi-tenant systems)

{' ' * (self.list_depth - 1)}- Processing untrusted user input or executing untrusted code
**Security risks:**

{' ' * (self.list_depth - 1)}- Agents can execute **arbitrary shell commands** with your user’s permissions

{' ' * (self.list_depth - 1)}- Agents can read any accessible file, including secrets (API keys, credentials, `.env` files)

{' ' * (self.list_depth - 1)}- Secrets may be exposed

{' ' * (self.list_depth - 1)}- File modifications and command execution are **permanent and irreversible**

{' ' * (self.list_depth - 1)}- Commands run directly on your host system

{' ' * (self.list_depth - 1)}- Commands can consume unlimited CPU, memory, disk
**Recommended safeguards:**

{' ' * (self.list_depth - 1)}- Enable [Human-in-the-Loop (HITL) middleware](/oss/python/deepagents/human-in-the-loop) to review and approve operations before execution. This is **strongly recommended**.

{' ' * (self.list_depth - 1)}- Run in dedicated development environments only. Never use on shared or production systems.

{' ' * (self.list_depth - 1)}- Use a [sandbox backend](/oss/python/deepagents/sandboxes) for production environments requiring shell execution.
**Note:** `virtual_mode=True` provides no security with shell access enabled, since commands can access any path on the system.

```
from deepagents import create_deep_agent
from deepagents.backends import LocalShellBackend

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=LocalShellBackend(root_dir=".", virtual_mode=True, env={"PATH": "/usr/bin:/bin"}),
)

```

**How it works:**

{' ' * (self.list_depth - 1)}- Extends `FilesystemBackend` with the `execute` tool for running shell commands on the host.

{' ' * (self.list_depth - 1)}- Commands run directly on your machine using `subprocess.run(shell=True)` with no sandboxing.

{' ' * (self.list_depth - 1)}- Supports `timeout` (default 120s), `max_output_bytes` (default 100,000), `env`, and `inherit_env` for environment variables.

{' ' * (self.list_depth - 1)}- Shell commands use `root_dir` as the working directory but can access any path on the system.

**Best for:**

{' ' * (self.list_depth - 1)}- Local coding assistants and development tools

{' ' * (self.list_depth - 1)}- Quick iteration during development when you trust the agent

### [​](#storebackend-langgraph-store)StoreBackend (LangGraph store)

```
from deepagents import create_deep_agent
from deepagents.backends import StoreBackend
from langgraph.store.memory import InMemoryStore

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=StoreBackend(
 namespace=lambda rt: (rt.server_info.user.identity,),
 ),
 store=InMemoryStore(), # Good for local dev; omit for LangSmith Deployment
)

```

When deploying to [LangSmith Deployment](/langsmith/deployment), omit the `store` parameter. The platform automatically provisions a store for your agent.
The `namespace` parameter controls data isolation. For multi-user deployments, always set a [namespace factory](/oss/python/deepagents/backends#namespace-factories) to isolate data per user or tenant.
**How it works:**

{' ' * (self.list_depth - 1)}- [`StoreBackend`](https://reference.langchain.com/python/deepagents/backends/store/StoreBackend) stores files in a LangGraph [`BaseStore`](https://reference.langchain.com/python/langchain-core/stores/BaseStore) provided by the runtime, enabling cross‑thread durable storage.

**Best for:**

{' ' * (self.list_depth - 1)}- When you already run with a configured LangGraph store (for example, Redis, Postgres, or cloud implementations behind [`BaseStore`](https://reference.langchain.com/python/langchain-core/stores/BaseStore)).

{' ' * (self.list_depth - 1)}- When you’re deploying your agent through [LangSmith Deployment](/langsmith/deployment) (a store is automatically provisioned for your agent).

#### [​](#namespace-factories)Namespace factories

A namespace factory controls where `StoreBackend` reads and writes data. It receives a LangGraph [`Runtime`](https://reference.langchain.com/python/langgraph/runtime/Runtime) and returns a tuple of strings used as the store namespace. Use namespace factories to isolate data between users, tenants, or assistants.
Pass the namespace factory to the `namespace` parameter when constructing a `StoreBackend`:

```
NamespaceFactory = Callable[[Runtime], tuple[str, ...]]

```

The `Runtime` provides:

{' ' * (self.list_depth - 1)}- `rt.context` — User-supplied context passed via LangGraph’s [context schema](https://langchain-ai.github.io/langgraph/concepts/runtime/) (for example, `user_id`)

{' ' * (self.list_depth - 1)}- `rt.server_info` — Server-specific metadata when running on LangGraph Server (assistant ID, graph ID, authenticated user)

{' ' * (self.list_depth - 1)}- `rt.execution_info` — Execution identity information (thread ID, run ID, checkpoint ID)

The `Runtime` argument is available in `deepagents>=0.5.2`. Earlier 0.5.x releases passed a `BackendContext` instead — see [migrating from `BackendContext`](#migrating-from-backendcontext) below. `rt.server_info` and `rt.execution_info` require `deepagents>=0.5.0`.
**Common namespace patterns:**

```
from deepagents.backends import StoreBackend

# Per-user: each user gets their own isolated storage
backend = StoreBackend(
 namespace=lambda rt: (rt.server_info.user.identity,),
)

# Per-assistant: all users of the same assistant share storage
backend = StoreBackend(
 namespace=lambda rt: (
 rt.server_info.assistant_id,
 ),
)

# Per-thread: storage scoped to a single conversation
backend = StoreBackend(
 namespace=lambda rt: (
 rt.execution_info.thread_id,
 ),
)

```

You can combine multiple components to create more specific scopes — for example, `(user_id, thread_id)` for per-user per-conversation isolation, or append a suffix like `"filesystem"` to disambiguate when the same scope uses multiple store namespaces.
Namespace components must contain only alphanumeric characters, hyphens, underscores, dots, `@`, `+`, colons, and tildes. Wildcards (`*`, `?`) are rejected to prevent glob injection.
The `namespace` parameter will be **required** in v0.5.0. Always set it explicitly for new code.
When no namespace factory is provided, the legacy default uses the `assistant_id` from LangGraph config metadata. This means all users of the same [assistant](/langsmith/assistants) share the same storage. For multi-user [going to production](/oss/python/deepagents/going-to-production), always provide a namespace factory.

### [​](#contexthubbackend)ContextHubBackend

```
from deepagents import create_deep_agent
from deepagents.backends import ContextHubBackend

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=ContextHubBackend("my-agent"),
)

```

`ContextHubBackend` stores files in a LangSmith Hub repo. Construct it with a repo identifier in `owner/name` or `name` format.
Set `LANGSMITH_API_KEY` before using `ContextHubBackend`.
**How it works:**

{' ' * (self.list_depth - 1)}- Pulls the Hub repo tree lazily on first use, then serves reads from an in-memory cache.

{' ' * (self.list_depth - 1)}- Persists writes and edits as Hub commits and updates the cache after successful commits.

{' ' * (self.list_depth - 1)}- Uses optimistic parent-commit writes (`parent_commit`): each push targets the latest known commit hash.

**Behavior and limits:**

{' ' * (self.list_depth - 1)}- If the repo does not exist, first pull is treated as empty; the first successful write can create the repo.

{' ' * (self.list_depth - 1)}- If another writer advances the repo first, your stale parent-commit write can fail. Re-pull and retry on conflict.

{' ' * (self.list_depth - 1)}- `upload_files()` accepts UTF-8 text. Non-UTF-8 files are rejected per path with `invalid_path`.

**Best for:**

{' ' * (self.list_depth - 1)}- LangSmith-native durable filesystem persistence without separately wiring a LangGraph `BaseStore`.

{' ' * (self.list_depth - 1)}- Workflows that benefit from Hub commit history on filesystem changes.

### [​](#compositebackend-router)CompositeBackend (router)

```
from deepagents import create_deep_agent
from deepagents.backends import CompositeBackend, StateBackend, StoreBackend
from langgraph.store.memory import InMemoryStore

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=CompositeBackend(
 default=StateBackend(),
 routes={
 "/memories/": StoreBackend(namespace=lambda _rt: ("memories",)),
 },
 ),
 store=InMemoryStore(), # Store passed to create_deep_agent, not backend
)

```

**How it works:**

{' ' * (self.list_depth - 1)}- [`CompositeBackend`](https://reference.langchain.com/python/deepagents/backends/composite/CompositeBackend) routes file operations to different backends based on path prefix.

{' ' * (self.list_depth - 1)}- Preserves the original path prefixes in listings and search results.

**Best for:**

{' ' * (self.list_depth - 1)}- When you want to give your agent both thread-scoped and cross-thread storage, a `CompositeBackend` allows you provide both a `StateBackend` and `StoreBackend`

{' ' * (self.list_depth - 1)}- When you have multiple sources of information that you want to provide to your agent as part of a single filesystem.

{' ' * (self.list_depth - 1)}- e.g. You have long-term memories stored under `/memories/` in one Store and you also have a custom backend that has documentation accessible at /docs/.

## [​](#specify-a-backend)Specify a backend

{' ' * (self.list_depth - 1)}- Pass a backend instance to `create_deep_agent(model=..., backend=...)`. The filesystem middleware uses it for all tooling.

{' ' * (self.list_depth - 1)}- The backend must implement `BackendProtocol` (for example, `StateBackend()`, `FilesystemBackend(root_dir=".")`, `StoreBackend()`, `ContextHubBackend("my-agent")`).

{' ' * (self.list_depth - 1)}- If omitted, the default is `StateBackend()`.

## [​](#route-to-different-backends)Route to different backends

Route parts of the namespace to different backends. Commonly used to persist `/memories/*` across threads and keep everything else thread-scoped.

```
from deepagents import create_deep_agent
from deepagents.backends import CompositeBackend, StateBackend, FilesystemBackend

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=CompositeBackend(
 default=StateBackend(),
 routes={
 "/memories/": FilesystemBackend(root_dir="/deepagents/myagent", virtual_mode=True),
 },
 )
)

```

Behavior:

{' ' * (self.list_depth - 1)}- `/workspace/plan.md` → `StateBackend` (thread-scoped)

{' ' * (self.list_depth - 1)}- `/memories/agent.md` → `FilesystemBackend` under `/deepagents/myagent`

{' ' * (self.list_depth - 1)}- `ls`, `glob`, `grep` aggregate results and show original path prefixes.

Notes:

{' ' * (self.list_depth - 1)}- Longer prefixes win (for example, route `"/memories/projects/"` can override `"/memories/"`).

{' ' * (self.list_depth - 1)}- For StoreBackend routing, ensure a store is provided via `create_deep_agent(model=..., store=...)` or provisioned by the platform.

{' ' * (self.list_depth - 1)}- Deep Agents write internal data (offloaded tool results, conversation history) to the default backend. Use `StateBackend` as the default to keep these artifacts ephemeral and avoid writing them to disk or a persistent store. See the [FilesystemBackend tip](#filesystembackend-local-disk) for a complete example.

## [​](#use-a-virtual-filesystem)Use a virtual filesystem

Build a custom backend to project a remote or database filesystem (e.g., S3 or Postgres) into the tools namespace.
Design guidelines:

{' ' * (self.list_depth - 1)}- 
Paths are absolute (`/x/y.txt`). Decide how to map them to your storage keys/rows.

{' ' * (self.list_depth - 1)}- 
Implement `ls` and `glob` efficiently (server-side filtering where available, otherwise local filter).

{' ' * (self.list_depth - 1)}- 
For external persistence (S3, Postgres, etc.), return `files_update=None` (Python) or omit `filesUpdate` (JS) in write/edit results — only in-memory state backends need to return a files update dict.

{' ' * (self.list_depth - 1)}- 
Use `ls` and `glob` as the method names.

{' ' * (self.list_depth - 1)}- 
Return structured result types with an `error` field for missing files or invalid patterns (do not raise).

S3-style outline:

```
from deepagents.backends.protocol import (
 BackendProtocol, WriteResult, EditResult, LsResult, ReadResult, GrepResult, GlobResult,
)

class S3Backend(BackendProtocol):
 def __init__(self, bucket: str, prefix: str = ""):
 self.bucket = bucket
 self.prefix = prefix.rstrip("/")

 def _key(self, path: str) -> str:
 return f"{self.prefix}{path}"

 def ls(self, path: str) -> LsResult:
 # List objects under _key(path); build FileInfo entries (path, size, modified_at)
 ...

 def read(self, file_path: str, offset: int = 0, limit: int = 2000) -> ReadResult:
 # Fetch object; return ReadResult(file_data=...) or ReadResult(error=...)
 ...

 def grep(self, pattern: str, path: str | None = None, glob: str | None = None) -> GrepResult:
 # Optionally filter server‑side; else list and scan content
 ...

 def glob(self, pattern: str, path: str = "/") -> GlobResult:
 # Apply glob relative to path across keys
 ...

 def write(self, file_path: str, content: str) -> WriteResult:
 # Enforce create‑only semantics; return WriteResult(path=file_path, files_update=None)
 ...

 def edit(self, file_path: str, old_string: str, new_string: str, replace_all: bool = False) -> EditResult:
 # Read → replace (respect uniqueness vs replace_all) → write → return occurrences
 ...

```

Postgres-style outline:

{' ' * (self.list_depth - 1)}- Table `files(path text primary key, content text, created_at timestamptz, modified_at timestamptz)`

{' ' * (self.list_depth - 1)}- Map tool operations onto SQL:

{' ' * (self.list_depth - 1)}- `ls` uses `WHERE path LIKE $1 || '%'`

{' ' * (self.list_depth - 1)}- `glob` filter in SQL or fetch then apply glob in Python

{' ' * (self.list_depth - 1)}- `grep` can fetch candidate rows by extension or last modified time, then scan lines

## [​](#permissions)Permissions

Use [permissions](/oss/python/deepagents/permissions) to declaratively control which files and directories the agent can read or write. Permissions apply to the built-in filesystem tools and are evaluated before the backend is called.

```
from deepagents import create_deep_agent, FilesystemPermission

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=CompositeBackend(
 default=StateBackend(),
 routes={
 "/memories/": StoreBackend(
 namespace=lambda rt: (rt.server_info.user.identity,),
 ),
 "/policies/": StoreBackend(
 namespace=lambda rt: (rt.context.org_id,),
 ),
 },
 ),
 permissions=[
 FilesystemPermission(
 operations=["write"],
 paths=["/policies/**"],
 mode="deny",
 ),
 ],
)

```

For the full set of options including rule ordering, subagent permissions, and composite backend interactions, see the [permissions guide](/oss/python/deepagents/permissions).

## [​](#add-policy-hooks)Add policy hooks

For custom validation logic beyond path-based allow/deny rules (rate limiting, audit logging, content inspection), enforce enterprise rules by subclassing or wrapping a backend.
Block writes/edits under selected prefixes (subclass):

```
from deepagents.backends.filesystem import FilesystemBackend
from deepagents.backends.protocol import WriteResult, EditResult

class GuardedBackend(FilesystemBackend):
 def __init__(self, *, deny_prefixes: list[str], **kwargs):
 super().__init__(**kwargs)
 self.deny_prefixes = [p if p.endswith("/") else p + "/" for p in deny_prefixes]

 def write(self, file_path: str, content: str) -> WriteResult:
 if any(file_path.startswith(p) for p in self.deny_prefixes):
 return WriteResult(error=f"Writes are not allowed under {file_path}")
 return super().write(file_path, content)

 def edit(self, file_path: str, old_string: str, new_string: str, replace_all: bool = False) -> EditResult:
 if any(file_path.startswith(p) for p in self.deny_prefixes):
 return EditResult(error=f"Edits are not allowed under {file_path}")
 return super().edit(file_path, old_string, new_string, replace_all)

```

Generic wrapper (works with any backend):

```
from deepagents.backends.protocol import (
 BackendProtocol, WriteResult, EditResult, LsResult, ReadResult, GrepResult, GlobResult,
)

class PolicyWrapper(BackendProtocol):
 def __init__(self, inner: BackendProtocol, deny_prefixes: list[str] | None = None):
 self.inner = inner
 self.deny_prefixes = [p if p.endswith("/") else p + "/" for p in (deny_prefixes or [])]

 def _deny(self, path: str) -> bool:
 return any(path.startswith(p) for p in self.deny_prefixes)

 def ls(self, path: str) -> LsResult:
 return self.inner.ls(path)

 def read(self, file_path: str, offset: int = 0, limit: int = 2000) -> ReadResult:
 return self.inner.read(file_path, offset=offset, limit=limit)
 def grep(self, pattern: str, path: str | None = None, glob: str | None = None) -> GrepResult:
 return self.inner.grep(pattern, path, glob)
 def glob(self, pattern: str, path: str = "/") -> GlobResult:
 return self.inner.glob(pattern, path)
 def write(self, file_path: str, content: str) -> WriteResult:
 if self._deny(file_path):
 return WriteResult(error=f"Writes are not allowed under {file_path}")
 return self.inner.write(file_path, content)
 def edit(self, file_path: str, old_string: str, new_string: str, replace_all: bool = False) -> EditResult:
 if self._deny(file_path):
 return EditResult(error=f"Edits are not allowed under {file_path}")
 return self.inner.edit(file_path, old_string, new_string, replace_all)

```

## [​](#migrate-from-backend-factories)Migrate from backend factories

The backend factory pattern is **deprecated** as of `deepagents` 0.5.0. Pass pre-constructed backend instances directly instead of factory functions.
Previously, backends like `StateBackend` and `StoreBackend` required a factory function that received a runtime object, because they needed runtime context (state, store) to operate. Backends now resolve this context internally via LangGraph’s `get_config()`, `get_store()`, and `get_runtime()` helpers, so you can pass instances directly.

### [​](#what-changed)What changed

| 
 | Before (deprecated) | After
 | `backend=lambda rt: StateBackend(rt)` | `backend=StateBackend()`
 | `backend=lambda rt: StoreBackend(rt)` | `backend=StoreBackend()`
 | `backend=lambda rt: CompositeBackend(default=StateBackend(rt), ...)` | `backend=CompositeBackend(default=StateBackend(), ...)`
 | `backend: (config) => new StateBackend(config)` | `backend: new StateBackend()`
 | `backend: (config) => new StoreBackend(config)` | `backend: new StoreBackend()`

### [​](#deprecated-apis)Deprecated APIs

| 
 | Deprecated | Replacement
 | Passing a callable to `backend=` in `create_deep_agent` | Pass a backend instance directly
 | `runtime` constructor argument on `StateBackend(runtime)` | `StateBackend()` (no arguments needed)
 | `runtime` constructor argument on `StoreBackend(runtime)` | `StoreBackend()` or `StoreBackend(namespace=..., store=...)`
 | `files_update` field on `WriteResult` and `EditResult` | State writes are now handled internally by the backend
 | `Command` wrapping in middleware write/edit tools | Tools return plain strings; no `Command(update=...)` needed
The factory pattern still works at runtime and emits a deprecation warning. Update your code to use direct instances before the next major version.

### [​](#migration-example)Migration example

```
# Before (deprecated)
from deepagents import create_deep_agent
from deepagents.backends import CompositeBackend, StateBackend, StoreBackend

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=lambda rt: CompositeBackend(
 default=StateBackend(rt),
 routes={"/memories/": StoreBackend(rt, namespace=lambda rt: (rt.server_info.user.identity,))},
 ),
)

# After
agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=CompositeBackend(
 default=StateBackend(),
 routes={"/memories/": StoreBackend(namespace=lambda rt: (rt.server_info.user.identity,))},
 ),
)

```

### [​](#migrating-from-backendcontext)Migrating from `BackendContext`

In `deepagents>=0.5.2` (Python) and `deepagents>=1.9.1` (TypeScript), namespace factories receive a LangGraph [`Runtime`](https://reference.langchain.com/python/langgraph/runtime/Runtime) directly instead of a `BackendContext` wrapper. The old `BackendContext` form still works via backwards-compatible `.runtime` and `.state` accessors, but those accessors emit a deprecation warning and will be removed in `deepagents>=0.7`.
**What changed:**

{' ' * (self.list_depth - 1)}- The factory argument is now a `Runtime`, not a `BackendContext`.

{' ' * (self.list_depth - 1)}- Drop the `.runtime` accessor — for example, `ctx.runtime.context.user_id` becomes `rt.server_info.user.identity`.

{' ' * (self.list_depth - 1)}- There is no direct replacement for `ctx.state`. Namespace info should be read-only and stable for the lifetime of a run, whereas state is mutable and changes step-to-step—deriving a namespace from it risks data ending up under inconsistent keys. If you have a use case that requires reading agent state, please [open an issue](https://github.com/langchain-ai/deepagents/issues).

```
# Before (deprecated, removed in v0.7)
StoreBackend(
 namespace=lambda ctx: (ctx.runtime.context.user_id,),
)

# After
StoreBackend(
 namespace=lambda rt: (rt.server_info.user.identity,),
)

```

## [​](#protocol-reference)Protocol reference

Backends must implement [`BackendProtocol`](https://reference.langchain.com/python/deepagents/backends/protocol/BackendProtocol).
Required methods:

{' ' * (self.list_depth - 1)}- `ls(path: str) -> LsResult`

{' ' * (self.list_depth - 1)}- Return entries with at least `path`. Include `is_dir`, `size`, `modified_at` when available. Sort by `path` for deterministic output.

{' ' * (self.list_depth - 1)}- `read(file_path: str, offset: int = 0, limit: int = 2000) -> ReadResult`

{' ' * (self.list_depth - 1)}- Return file data on success. On missing file, return `ReadResult(error="Error: File '/x' not found")`.

{' ' * (self.list_depth - 1)}- `grep(pattern: str, path: Optional[str] = None, glob: Optional[str] = None) -> GrepResult`

{' ' * (self.list_depth - 1)}- Return structured matches. On error, return `GrepResult(error="...")` (do not raise).

{' ' * (self.list_depth - 1)}- `glob(pattern: str, path: str = "/") -> GlobResult`

{' ' * (self.list_depth - 1)}- Return matched files as `FileInfo` entries (empty list if none).

{' ' * (self.list_depth - 1)}- `write(file_path: str, content: str) -> WriteResult`

{' ' * (self.list_depth - 1)}- Create-only. On conflict, return `WriteResult(error=...)`. On success, set `path` and for state backends set `files_update={...}`; external backends should use `files_update=None`.

{' ' * (self.list_depth - 1)}- `edit(file_path: str, old_string: str, new_string: str, replace_all: bool = False) -> EditResult`

{' ' * (self.list_depth - 1)}- Enforce uniqueness of `old_string` unless `replace_all=True`. If not found, return error. Include `occurrences` on success.

Supporting types:

{' ' * (self.list_depth - 1)}- `LsResult(error, entries)` — `entries` is a `list[FileInfo]` on success, `None` on failure.

{' ' * (self.list_depth - 1)}- `ReadResult(error, file_data)` — `file_data` is a `FileData` dict on success, `None` on failure.

{' ' * (self.list_depth - 1)}- `GrepResult(error, matches)` — `matches` is a `list[GrepMatch]` on success, `None` on failure.

{' ' * (self.list_depth - 1)}- `GlobResult(error, matches)` — `matches` is a `list[FileInfo]` on success, `None` on failure.

{' ' * (self.list_depth - 1)}- `WriteResult(error, path, files_update)`

{' ' * (self.list_depth - 1)}- `EditResult(error, path, files_update, occurrences)`

{' ' * (self.list_depth - 1)}- `FileInfo` with fields: `path` (required), optionally `is_dir`, `size`, `modified_at`.

{' ' * (self.list_depth - 1)}- `GrepMatch` with fields: `path`, `line`, `text`.

{' ' * (self.list_depth - 1)}- `FileData` with fields: `content` (str), `encoding` (`"utf-8"` or `"base64"`), `created_at`, `modified_at`.
:::

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/backends.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Context engineering in Deep AgentsPrevious](/oss/python/deepagents/context-engineering)[SubagentsNext](/oss/python/deepagents/subagents)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
