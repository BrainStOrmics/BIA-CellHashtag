# Permissions

[Core capabilities](/oss/python/deepagents/harness)

# Permissions
Copy page

Control filesystem access with declarative permission rules for Deep AgentsCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Control which files and directories an agent can read or write to using declarative permission rules. Pass a list of rules to `permissions=` and the agent’s built-in filesystem tools respect them.
Permissions require `deepagents>=0.5.2`.
Permissions only apply to the built-in filesystem tools (`ls`, `read_file`, `glob`, `grep`, `write_file`, `edit_file`). Custom tools and MCP tools that access the filesystem are not covered. Permissions also do not apply to [sandbox backends](/oss/python/deepagents/sandboxes), which support arbitrary command execution via the `execute` tool.
Use `permissions` when you need **path-based allow/deny rules** on the built-in filesystem tools. Use [backend policy hooks](/oss/python/deepagents/backends#add-policy-hooks) when you need custom validation logic (rate limiting, audit logging, content inspection) or need to control custom tools.

## [​](#basic-usage)Basic usage

Pass a list of [`FilesystemPermission`](https://reference.langchain.com/python/deepagents/middleware/permissions/FilesystemPermission) rules to [`create_deep_agent`](https://reference.langchain.com/python/deepagents/graph/create_deep_agent). Rules are evaluated in declaration order. The first matching rule wins. If no rule matches, the operation is allowed.

```
from deepagents import FilesystemPermission, create_deep_agent

# Read-only agent: deny all writes
agent = create_deep_agent(
 model=model,
 backend=backend,
 permissions=[
 FilesystemPermission(
 operations=["write"],
 paths=["/**"],
 mode="deny",
 ),
 ],
)

```

## [​](#rule-structure)Rule structure

Each `FilesystemPermission` has three fields:

| 
 | Field | Type | Description
 | `operations` | `list["read" | "write"]` | Operations this rule applies to. `"read"` covers `ls`, `read_file`, `glob`, `grep`. `"write"` covers `write_file`, `edit_file`.
 | `paths` | `list[str]` | Glob patterns for matching file paths (e.g., `["/workspace/**"]`). Supports `**` for recursive matching and `{a,b}` for alternation.
 | `mode` | `"allow" | "deny"` | Whether to allow or deny matching operations. Defaults to `"allow"`.
Rules use first-match-wins evaluation: the first rule whose `operations` and `paths` match the current call determines the outcome. If no rule matches, the call is **allowed** (permissive default).

## [​](#examples)Examples

### [​](#isolate-to-a-workspace-directory)Isolate to a workspace directory

Allow reads and writes only under `/workspace/` and deny everything else:

```
agent = create_deep_agent(
 model=model,
 backend=backend,
 permissions=[
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/workspace/**"],
 mode="allow",
 ),
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/**"],
 mode="deny",
 ),
 ],
)

```

### [​](#protect-specific-files)Protect specific files

```
agent = create_deep_agent(
 model=model,
 backend=backend,
 permissions=[
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/workspace/.env", "/workspace/examples/**"],
 mode="deny",
 ),
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/workspace/**"],
 mode="allow",
 ),
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/**"],
 mode="deny",
 ),
 ],
)

```

### [​](#read-only-memory)Read-only memory

Allow the agent to read memory files but prevent it from modifying them. This is useful for organization-wide policies or shared knowledge bases that should only be updated by application code. See [read-only vs writable memory](/oss/python/deepagents/memory#read-only-vs-writable-memory) for more context.

```
from deepagents.backends import CompositeBackend, StateBackend, StoreBackend

agent = create_deep_agent(
 model=model,
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
 paths=["/memories/**", "/policies/**"],
 mode="deny",
 ),
 ],
)

```

### [​](#deny-all-access)Deny all access

Block all reads and writes. This is a restrictive baseline you can layer more specific allow rules on top of:

```
agent = create_deep_agent(
 model=model,
 backend=backend,
 permissions=[
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/**"],
 mode="deny",
 ),
 ],
)

```

### [​](#rule-ordering)Rule ordering

Because of first-match-wins, rule order matters. Place more specific rules before broader ones:

```
# Correct: deny .env, allow workspace, deny everything else
correct_permissions = [
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/workspace/.env"],
 mode="deny",
 ),
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/workspace/**"],
 mode="allow",
 ),
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/**"],
 mode="deny",
 ),
]

# Bug: /workspace/** matches .env first, so the deny never triggers
incorrect_permissions = [
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/workspace/**"],
 mode="allow",
 ),
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/workspace/.env"],
 mode="deny", # never reached
 ),
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/**"],
 mode="deny",
 ),
]

```

## [​](#subagent-permissions)Subagent permissions

[Subagents](/oss/python/deepagents/subagents) inherit the parent agent’s permissions by default. To give a subagent different permissions, set the `permissions` field in its spec. This **replaces** the parent’s rules entirely.

```
agent = create_deep_agent(
 model=model,
 backend=backend,
 permissions=[
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/workspace/**"],
 mode="allow",
 ),
 FilesystemPermission(
 operations=["read", "write"],
 paths=["/**"],
 mode="deny",
 ),
 ],
 subagents=[
 {
 "name": "auditor",
 "description": "Read-only code reviewer",
 "system_prompt": "Review the code for issues.",
 "permissions": [
 FilesystemPermission(
 operations=["write"],
 paths=["/**"],
 mode="deny",
 ),
 FilesystemPermission(
 operations=["read"],
 paths=["/workspace/**"],
 mode="allow",
 ),
 FilesystemPermission(
 operations=["read"],
 paths=["/**"],
 mode="deny",
 ),
 ],
 }
 ],
)

```

## [​](#composite-backends)Composite backends

When using a [`CompositeBackend`](https://reference.langchain.com/python/deepagents/backends/composite/CompositeBackend) with a sandbox default, every permission path must be scoped under a known route prefix. Sandboxes support arbitrary command execution, so path-based restrictions alone cannot prevent filesystem access through shell commands. Scoping permissions to route-specific [backends](/oss/python/deepagents/backends) avoids this conflict.

```
from deepagents.backends import CompositeBackend

composite = CompositeBackend(
 default=sandbox,
 routes={"/memories/": memories_backend},
)

# Works: permissions are scoped to the /memories/ route
agent = create_deep_agent(
 model=model,
 backend=composite,
 permissions=[
 FilesystemPermission(
 operations=["write"],
 paths=["/memories/**"],
 mode="deny",
 ),
 ],
)

```

Permissions that include paths outside any route raise `NotImplementedError`:

```
# Raises NotImplementedError: /workspace/** hits the sandbox default
try:
 create_deep_agent(
 model=model,
 backend=composite,
 permissions=[
 FilesystemPermission(
 operations=["write"],
 paths=["/workspace/**"],
 mode="deny",
 ),
 ],
 )
except NotImplementedError:
 pass

# Also raises: /** covers both routes and the default
try:
 create_deep_agent(
 model=model,
 backend=composite,
 permissions=[
 FilesystemPermission(
 operations=["read"],
 paths=["/**"],
 mode="deny",
 ),
 ],
 )
except NotImplementedError:
 pass

```

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/permissions.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Human-in-the-loopPrevious](/oss/python/deepagents/human-in-the-loop)[MemoryNext](/oss/python/deepagents/memory)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
