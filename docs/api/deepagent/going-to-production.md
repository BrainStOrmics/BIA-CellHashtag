# Going to production

[Deployment](/langsmith/deploy-managed-deep-agent)

# Going to production
Copy page

Take your deep agent to production with persistent memory, sandboxes, resilience middleware, and deployment optionsCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.This guide covers considerations for taking a deep agent from a local prototype to a production deployment. It walks through scoping memory, configuring execution environments, adding guardrails, and connecting a frontend.

## [​](#overview)Overview

Agents use information from memory and their execution environment to accomplish tasks.
In production, there are a few primitives that determine how information is shared and accessed:

{' ' * (self.list_depth - 1)}- **Thread**: a single conversation. Message history and scratch files are scoped to the thread by default and don’t carry over.

{' ' * (self.list_depth - 1)}- **User**: someone interacting with your agent. Memory and files can be private to a user or shared across users. Identity and authorization comes from your [auth layer](/langsmith/auth).

{' ' * (self.list_depth - 1)}- **Assistant**: a configured agent instance. Memory and files can be tied to one assistant or shared across all of them.

This page covers:

{' ' * (self.list_depth - 1)}- **[LangSmith Deployments](#langsmith-deployments)**: managed infrastructure with auth, webhooks, and cron

{' ' * (self.list_depth - 1)}- **[Production considerations](#production-considerations)**: invocation, multi-tenancy, authentication, credentials, async, and durability

{' ' * (self.list_depth - 1)}- **[Memory](#memory)**: persist information across conversations

{' ' * (self.list_depth - 1)}- **[Execution environment](#execution-environment)**: file storage and code execution

{' ' * (self.list_depth - 1)}- **[Guardrails](#guardrails)**: rate limiting, error handling, and data privacy

{' ' * (self.list_depth - 1)}- **[Frontend](#frontend)**: connect your UI to a deployed agent

## [​](#langsmith-deployments)LangSmith Deployments

![Managed Deep Agents packages your agent configuration, tools, and runtime settings for LangSmith](https://mintcdn.com/langchain-5e9cc07a/9oyV6nbtSbBRfaE1/oss/images/deepagents/production/deepagents-deploy-config.png?fit=max&auto=format&n=9oyV6nbtSbBRfaE1&q=85&s=5c3f2961994afe5fe67a2f5c9e9ba7ac)
The recommended path for taking a Deep Agent to production is [Managed Deep Agents](/langsmith/deploy-managed-deep-agent), an API-first hosted runtime for creating, running, and operating deep agents in LangSmith. Managed Deep Agents is currently in private preview ([join the waitlist](https://www.langchain.com/langsmith-managed-deep-agents-waitlist)). For teams that need custom application code, custom routes, advanced authentication, or full Agent Server APIs, you can configure a [LangSmith Deployment](/langsmith/deployment) directly. Either path provisions the infrastructure your agent needs: [threads](/langsmith/use-threads), [runs](/langsmith/runs), a store, and a checkpointer, so you don’t have to set these up yourself. A traditional LangSmith Deployment also gives you [authentication](/langsmith/auth), [webhooks](/langsmith/use-webhooks), [cron jobs](/langsmith/cron-jobs), and [observability](/langsmith/observability) out of the box, and can expose your agent via [MCP](/langsmith/server-mcp) or [A2A](/langsmith/server-a2a).
All code snippets on this page use the following `langgraph.json` unless otherwise specified:
langgraph.json
```
{
 "dependencies": ["."],
 "graphs": {
 "agent": "./agent.py:agent"
 },
 "env": ".env"
}

```

`langgraph.json` is the configuration file that tells the LangGraph platform how to build and run your application. It lives at the root of your project and is required for both local development (with `langgraph dev`) and production deployment. The key fields are:

| 
 | Field | Description
 | `dependencies` | Packages to install. `["."]` installs the current directory as a package (reads from `requirements.txt`, `pyproject.toml`, or `package.json`).
 | `graphs` | Maps graph IDs to their code locations. Each entry is `"<id>": "./<file>:<variable>"`, where `<id>` is the name you use to invoke the graph via the API, and `<variable>` is the compiled graph or constructor function exported from `<file>`.
 | `env` | Path to a `.env` file with environment variables (API keys, secrets). These are set at build time and available at runtime.
For the full set of configuration options (custom Docker steps, store indexing, auth handlers, and more), see [application structure](/oss/python/langgraph/application-structure).

## [​](#production-considerations)Production considerations

### [​](#invoking-the-agent)Invoking the agent

In production, every invocation should carry two run-level parameters:

{' ' * (self.list_depth - 1)}- **`thread_id`** (passed via `config={"configurable": {"thread_id": ...}}`): a stable identifier for the conversation. The [checkpointer](#durability) uses it to persist and resume message history, so follow-up turns continue the same conversation. Generate a new `thread_id` to start a fresh conversation.

{' ' * (self.list_depth - 1)}- **`context`**: per-run data your tools and middleware read at invocation time, for example `user_id`, API keys, feature flags, or session metadata. Define the shape with `context_schema` and access it via `runtime.context`. See [Runtime context](/oss/python/deepagents/context-engineering#runtime-context).

The two are independent and almost always passed together:
GoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from dataclasses import dataclass

from deepagents import create_deep_agent
from langchain_core.utils.uuid import uuid7

@dataclass
class Context:
 user_id: str

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 context_schema=Context,
)

# Start a conversation
config = {"configurable": {"thread_id": str(uuid7())}}
agent.invoke(
 {"messages": [{"role": "user", "content": "Plan a 3-day trip to Tokyo"}]},
 config=config,
 context=Context(user_id="user-123"),
)

# Follow-up on the same conversation: reuse the same thread_id
agent.invoke(
 {"messages": [{"role": "user", "content": "Make it 5 days instead"}]},
 config=config,
 context=Context(user_id="user-123"),
)

```

When deploying with the LangGraph SDK, the SDK manages threads for you and you pass the returned `thread_id` to each run:

```
from langgraph_sdk import get_client

client = get_client(url="<DEPLOYMENT_URL>", api_key="<LANGSMITH_API_KEY>")

thread = await client.threads.create()
async for chunk in client.runs.stream(
 thread["thread_id"],
 "agent",
 input={"messages": [{"role": "user", "content": "Plan a 3-day trip to Tokyo"}]},
 context={"user_id": "user-123"},
 stream_mode="updates",
):
 print(chunk.data)

```

`thread_id` scopes the *conversation* (message history, checkpoints). `context` carries *per-run* data your tools and middleware read. They are independent: changing one does not affect the other, and you can pass either or both.

### [​](#multi-tenancy)Multi-tenancy

When your agent serves multiple users, you need to handle three concerns: verifying who each user is, controlling what they can access, and managing the credentials the agent uses to act on their behalf.

![Three authentication layers compose: end-user auth, agent-acting-as-user auth, and team RBAC](https://mintcdn.com/langchain-5e9cc07a/9oyV6nbtSbBRfaE1/oss/images/deepagents/production/auth-layers.png?fit=max&auto=format&n=9oyV6nbtSbBRfaE1&q=85&s=11137d8bae74093c459416fd3e541d68)

#### [​](#user-identity-and-access-control)User identity and access control

[LangSmith Deployments](/langsmith/deployment) supports [custom authentication](/langsmith/custom-auth) to establish user identity and [authorization handlers](/langsmith/auth) to control access to resources like threads, assistants, and store namespaces. Authorization handlers run after authentication succeeds and can:

{' ' * (self.list_depth - 1)}- Tag resources with ownership metadata (e.g., `owner: user_id`)

{' ' * (self.list_depth - 1)}- Return filters so users only see their own resources

{' ' * (self.list_depth - 1)}- Deny access with HTTP 403 for unauthorized operations

For a step-by-step tutorial, see [Make conversations private](/langsmith/resource-auth). For a walkthrough, watch the [custom auth video](https://www.youtube.com/watch?v=DkNqgCz8cjE).
How you [scope memory](#scoping) and [execution environments](#execution-environment) determines what data is shared between users. See the sections below for details.

#### [​](#team-access-control-rbac)Team access control (RBAC)

LangSmith’s [role-based access control](/langsmith/rbac) governs who on your team can deploy, configure, and monitor agents. This is separate from end-user authorization above.

| 
 | Role | Access
 | Workspace Admin | Full permissions including settings and member management
 | Workspace Editor | Create and modify resources, but cannot delete runs or manage members
 | Workspace Viewer | Read-only access
Custom roles with granular permissions are available on Enterprise plans. See the [RBAC reference](/langsmith/rbac) for the full permission model.

#### [​](#end-user-credentials)End-user credentials

When your agent needs to call external APIs on behalf of a user (e.g., reading their GitHub repos, sending Slack messages, querying their data warehouse), you need a way to pass the user’s credentials through to the agent without hardcoding them.
**OAuth via Agent Auth.** [Agent Auth](/langsmith/agent-auth) provides a managed OAuth 2.0 flow. Configure an OAuth provider, and the agent can request tokens scoped to each user. On first use, the agent [interrupts](/oss/python/langgraph/interrupts) execution and presents an OAuth consent URL. After the user authenticates, the agent resumes with a valid token. Tokens are stored and refreshed automatically.

```
from langchain_auth import Client
from langchain.tools import tool, ToolRuntime

auth_client = Client()

# Inside your agent's tool:
@tool
async def github_action(runtime: ToolRuntime):
 """Perform an action on behalf of the user via GitHub."""
 auth_result = await auth_client.authenticate(
 provider="github",
 scopes=["repo", "read:org"],
 user_id=runtime.server_info.user.identity,
 )
 # Use auth_result.token for GitHub API calls on the user's behalf

```

**Credential injection for sandboxes.** If your agent runs code inside a [sandbox](#sandboxes) that calls external APIs, the [sandbox auth proxy](/langsmith/sandbox-auth-proxy) can inject credentials into outbound requests automatically, so sandbox code never receives raw API keys. See [Managing secrets](#managing-secrets) for setup details.
**Workspace secrets.** For API keys shared across all users (for example your organization’s LLM provider keys, search API keys), store them as [workspace secrets](/langsmith/set-up-hierarchy#configure-workspace-settings) in LangSmith. See [Managing secrets](#managing-secrets) for details.

### [​](#async)Async

LLM-based applications are heavily I/O-bound: calling language models, databases, and external services. Async programming lets these operations run concurrently instead of blocking, improving throughput and responsiveness.
LangChain follows the convention of prefixing `a` to async method names (e.g., `ainvoke`, `abefore_agent`, `astream`). Sync and async variants live in the same class or namespace.
When building for production:

{' ' * (self.list_depth - 1)}- **Create async tools.** LangChain runs sync tools in a separate thread to avoid blocking, but native async avoids the threading overhead entirely.

{' ' * (self.list_depth - 1)}- **Use async middleware methods.** Custom [middleware](/oss/python/langchain/middleware/custom) should implement async hooks (e.g., `abefore_agent` instead of `before_agent`).

{' ' * (self.list_depth - 1)}- **Use async for external resource lifecycle.** Creating [sandboxes](#sandboxes) or connecting to [MCP servers](/oss/python/langchain/mcp) involves network calls and should be awaited. This is why [graph factories](/langsmith/graph-rebuild) that provision these resources are async.

### [​](#durability)Durability

Deep Agents run on LangGraph, which provides [durable execution](/oss/python/langgraph/durable-execution) out of the box. The [persistence](/oss/python/langgraph/persistence) layer checkpoints state at each step, so a run interrupted by a failure, timeout, or [human-in-the-loop](/oss/python/langgraph/interrupts) pause resumes from its last recorded state without reprocessing previous steps. For long-running deep agents that spawn many subagents, this means a mid-run failure doesn’t lose completed work.

![Durable execution: when a worker crashes mid-run, another worker picks the run up from the latest checkpoint](https://mintcdn.com/langchain-5e9cc07a/9oyV6nbtSbBRfaE1/oss/images/deepagents/production/durable-execution.png?fit=max&auto=format&n=9oyV6nbtSbBRfaE1&q=85&s=d3bfd69460769dba142c68c7a20ae43b)
Checkpointing also enables:

{' ' * (self.list_depth - 1)}- **Indefinite [interrupts](/oss/python/langgraph/interrupts).** Human-in-the-loop workflows can pause for minutes or days and resume exactly where they left off.

{' ' * (self.list_depth - 1)}- **[Time travel](/oss/python/langgraph/use-time-travel).** Every checkpointed step is a snapshot you can rewind to, letting you replay from an earlier state if something goes wrong.

{' ' * (self.list_depth - 1)}- **Safe handling of sensitive operations.** For workflows involving payments or other irreversible actions, checkpoints provide an audit trail and a recovery point to inspect the exact state that led to an action.

[LangSmith Deployments](/langsmith/deployment) configure a persistent checkpointer automatically. If you are self-hosting, see [persistence](/oss/python/langgraph/persistence) for setup instructions.

## [​](#memory)Memory

Without memory, every conversation starts from scratch. Memory lets your agent retain information across conversations (user preferences, learned instructions, past experiences) so it can personalize its behavior over time. For an overview of memory types, see the [memory concepts guide](/oss/python/concepts/memory).

![Short-term memory is scoped to a single thread via checkpoints; long-term memory persists across threads via the store](https://mintcdn.com/langchain-5e9cc07a/9oyV6nbtSbBRfaE1/oss/images/deepagents/production/memory.png?fit=max&auto=format&n=9oyV6nbtSbBRfaE1&q=85&s=3ec585271dcd8d62e0207d79d68c296b)

### [​](#scoping)Scoping

Memory is always persistent across conversations. The main question is how it’s scoped across user and assistant boundaries. The right scope depends on who should see and modify the data:

| 
 | Scope | Namespace | Use case | Example
 | **User** (recommended default) | `(user_id)` | Per-user preferences and context | ”I prefer concise responses”
 | **Assistant** | `(assistant_id)` | Shared instructions for one assistant | ”Cap posts at 280 characters”
 | **Global** | `(org_id)` | Read-only policies for all users and assistants | ”Never disclose internal pricing”
Shared memory (assistant, user, or organization scope) is a vector for prompt injection. If one user can write to memory that another user’s conversation reads, a malicious user could inject instructions into that shared state. Enforce read-only access where appropriate. For example, make organization-wide policies writable only through application code, not by the agent itself. Use [permissions](/oss/python/deepagents/permissions) to declaratively deny writes to shared paths, or [backend policy hooks](/oss/python/deepagents/backends#add-policy-hooks) for custom validation logic.

### [​](#configuration)Configuration

In Deep Agents, memory is stored as files in a virtual filesystem. By default, files are scoped to a single thread (conversation) and not shared across threads.
Otherwise, to share memory across threads, route a path like `/memories/` to a [StoreBackend](https://reference.langchain.com/python/deepagents/backends/store/StoreBackend) that writes to the LangGraph [Store](/langsmith/custom-store). Use a [CompositeBackend](https://reference.langchain.com/python/deepagents/backends/composite/CompositeBackend) to give the agent both thread-scoped scratch space and cross-thread [long-term memory](/oss/python/deepagents/memory).
The `rt.server_info` and `rt.execution_info` namespace patterns shown below require `deepagents>=0.5.0`.

{' ' * (self.list_depth - 1)}- User (recommended)
{' ' * (self.list_depth - 1)}- Assistant
{' ' * (self.list_depth - 1)}- User
{' ' * (self.list_depth - 1)}- OrganizationNamespace by `user_id`. Each user gets their own private memory. This is the recommended default since most applications deploy a single assistant.agent.py
```
from deepagents import create_deep_agent
from deepagents.backends import CompositeBackend, StateBackend, StoreBackend

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=CompositeBackend(
 default=StateBackend(),
 routes={
 "/memories/": StoreBackend(
 namespace=lambda rt: (
 rt.server_info.assistant_id,
 rt.server_info.user.identity,
 ),
 ),
 },
 ),
 system_prompt="""You have persistent memory at /memories/.

 Read /memories/instructions.txt at the start of each conversation for
 accumulated knowledge and preferences. When you learn something that
 should persist, update that file.""",
)

```
Namespace by `assistant_id`. Memory is shared across all users of the same assistant, so any user can read or update it. Use this for shared instructions or knowledge that applies to everyone using a given assistant (e.g., “always reply in formal tone”).agent.py
```
from deepagents import create_deep_agent
from deepagents.backends import CompositeBackend, StateBackend, StoreBackend

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=CompositeBackend(
 default=StateBackend(),
 routes={
 "/memories/": StoreBackend(
 namespace=lambda rt: (
 rt.server_info.assistant_id,
 ),
 ),
 },
 ),
)

```
Namespace by `user_id` alone. Memory follows the user across all assistants. Use this for a global user profile (name, timezone, communication preferences) that should apply regardless of which assistant the user is talking to.agent.py
```
from deepagents import create_deep_agent
from deepagents.backends import CompositeBackend, StateBackend, StoreBackend

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=CompositeBackend(
 default=StateBackend(),
 routes={
 "/memories/": StoreBackend(
 namespace=lambda rt: (rt.server_info.user.identity,),
 ),
 },
 ),
)

```
Namespace by `org_id`. Memory is shared across all users and all assistants. Typically used for organization-wide policies (compliance rules, brand guidelines) that should be read-only for the agent. Write access should be restricted to application code to prevent prompt injection.agent.py
```
from deepagents import create_deep_agent
from deepagents.backends import CompositeBackend, StateBackend, StoreBackend

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=CompositeBackend(
 default=StateBackend(),
 routes={
 "/memories/": StoreBackend(
 namespace=lambda rt: (rt.context.org_id,),
 ),
 },
 ),
)

```

You can also read and write to the store from your application code using the [Store API](/langsmith/custom-store). See [Advanced usage](/oss/python/deepagents/memory#advanced-usage) for examples.
For the full namespace factory API, see [namespace factories](/oss/python/deepagents/backends#namespace-factories). For memory patterns like self-improving instructions and knowledge bases, see [long-term memory](/oss/python/deepagents/memory).

## [​](#execution-environment)Execution environment

Locally, agents can read and write files on disk and run shell commands directly. In production, you need to think about isolation and persistence. The right setup depends on whether your agent needs to execute code:

{' ' * (self.list_depth - 1)}- **Filesystem backends** are enough if your agent only reads and writes files. Choose a backend that matches your persistence needs: thread-scoped scratch space, cross-thread storage, or a mix of both.

{' ' * (self.list_depth - 1)}- **Sandboxes** add an isolated container with an `execute` tool for running shell commands. Use a sandbox if your agent needs to run code, install packages, or do anything beyond file I/O.

### [​](#filesystem)Filesystem

Choose a backend based on what needs to persist:

{' ' * (self.list_depth - 1)}- 
[StateBackend](https://reference.langchain.com/python/deepagents/backends/state/StateBackend) (default): thread-scoped scratch space. Files persist across turns within a thread via your checkpointer but are not shared across threads. Checkpointed at every step, so avoid writing large files.

{' ' * (self.list_depth - 1)}- 
[StoreBackend](https://reference.langchain.com/python/deepagents/backends/store/StoreBackend): cross-thread storage that survives across conversations. Scope with a [namespace factory](/oss/python/deepagents/backends#namespace-factories).

{' ' * (self.list_depth - 1)}- 
[CompositeBackend](https://reference.langchain.com/python/deepagents/backends/composite/CompositeBackend): mix both. Thread-scoped scratch space by default with cross-thread routes for specific paths like `/memories/`.

{' ' * (self.list_depth - 1)}- 
[`ContextHubBackend`](/oss/python/deepagents/backends#contexthubbackend): durable files in a LangSmith Hub repo (`owner/name` or `name`). Use this when you want LangSmith-native persistence without provisioning a separate LangGraph store.

For the full list of backends and how to build custom ones, see [backends](/oss/python/deepagents/backends).
`FilesystemBackend` and `LocalShellBackend` access the host directly. Don’t use them in deployed agents.

### [​](#sandboxes)Sandboxes

If your agent needs to run code (not just read and write files), use a [sandbox](/oss/python/deepagents/sandboxes). Sandboxes provide both a filesystem and an `execute` tool for running shell commands, all inside an isolated container. This isolation also protects your host: if the agent’s code exhausts memory or crashes, only the sandbox is affected. Your server keeps running.

#### [​](#lifecycle)Lifecycle

The key decision is how long a sandbox lives. Does each conversation get a fresh one, or do conversations share a persistent environment?

| 
 | Scope | Sandbox ID stored on | Lifecycle | Example use case
 | **Thread-scoped** | [Thread](/langsmith/use-threads) metadata | Fresh per conversation, cleaned up on TTL | A data analysis bot where each conversation starts clean
 | **Assistant-scoped** | [Assistant](/langsmith/assistants) config | Shared across all conversations | A coding assistant that maintains a cloned repo across conversations
The examples below use an async [graph factory](/langsmith/graph-rebuild) instead of a static graph because the sandbox needs the `thread_id` or `assistant_id` to look up or create the correct sandbox. Graph factories don’t receive a full `Runtime` (no `server_info` or `execution_info`); instead, accept a `RunnableConfig` and read `thread_id` and `assistant_id` from `config["configurable"]`. The factory is async because sandbox creation is an I/O-bound operation that requires per-run information only available at invocation time.

{' ' * (self.list_depth - 1)}- Thread-scoped (most common)
{' ' * (self.list_depth - 1)}- Assistant-scopedEach conversation gets its own sandbox. The [graph factory](/langsmith/graph-rebuild) reads `thread_id` from the run config, so each [thread](/langsmith/use-threads) automatically gets its own isolated environment. The provider’s label-based lookup handles deduplication across runs. Cleaned up when the sandbox [TTL](/langsmith/configure-ttl) expires.agent.py
```
from daytona import CreateSandboxFromSnapshotParams, Daytona
from deepagents import create_deep_agent
from langchain_core.runnables import RunnableConfig
from langchain_daytona import DaytonaSandbox

client = Daytona()

async def agent(config: RunnableConfig):
 thread_id = config["configurable"]["thread_id"]
 try:
 sandbox = await client.find_one(labels={"thread_id": thread_id})
 except Exception:
 sandbox = await client.create(
 CreateSandboxFromSnapshotParams(
 labels={"thread_id": thread_id},
 auto_delete_interval=3600, # TTL: clean up when idle
 )
 )
 return create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=DaytonaSandbox(sandbox=sandbox)
 )

```
All conversations share one sandbox. The [graph factory](/langsmith/graph-rebuild) reads the [assistant](/langsmith/assistants) ID from `config["configurable"]`, so every thread on the same assistant returns to the same environment. Files, installed packages, and cloned repositories persist across conversations.agent.py
```
from daytona import CreateSandboxFromSnapshotParams, Daytona
from deepagents import create_deep_agent
from langchain_core.runnables import RunnableConfig
from langchain_daytona import DaytonaSandbox

client = Daytona()

async def agent(config: RunnableConfig):
 assistant_id = config["configurable"]["assistant_id"]
 try:
 sandbox = await client.find_one(labels={"assistant_id": assistant_id})
 except Exception:
 sandbox = await client.create(
 CreateSandboxFromSnapshotParams(labels={"assistant_id": assistant_id})
 )
 return create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=DaytonaSandbox(sandbox=sandbox)
 )

```
Assistant-scoped sandboxes accumulate files, installed packages, and other in-sandbox state over time. Configure a TTL with your sandbox provider, use snapshots to reset periodically, or implement cleanup logic to prevent the sandbox’s disk and memory from growing unbounded.
Because the `agent` variable is an async function (not a compiled graph), the server treats it as a [graph factory](/langsmith/graph-rebuild) and calls it on each run, injecting the config. The factory looks up or creates the sandbox via the provider’s label-based search and returns a fresh agent graph wired to that sandbox.
Once deployed with `langgraph deploy`, invoke the agent from your application code using the SDK. The client-side code is the same regardless of scope. The scoping is handled entirely in the agent factory above, but the behavior differs:

{' ' * (self.list_depth - 1)}- Thread-scoped
{' ' * (self.list_depth - 1)}- Assistant-scopedEach thread gets its own sandbox. Follow-up messages within the same thread reuse the same sandbox, but a new thread always starts fresh with no leftover files or installed packages from previous conversations.client.py
```
from langgraph_sdk import get_client

client = get_client(url="<DEPLOYMENT_URL>", api_key="<LANGSMITH_API_KEY>")

# Conversation 1: install pandas and analyze data
thread_1 = await client.threads.create()
async for chunk in client.runs.stream(
 thread_1["thread_id"],
 "agent",
 input={"messages": [{"role": "human", "content": "Install pandas and analyze sales_data.csv"}]},
 stream_mode="updates",
):
 print(chunk.data)

# Follow-up in the same conversation — pandas is still installed
async for chunk in client.runs.stream(
 thread_1["thread_id"],
 "agent",
 input={"messages": [{"role": "human", "content": "Now plot the results"}]},
 stream_mode="updates",
):
 print(chunk.data)

# Conversation 2: fresh sandbox — pandas is NOT installed, no files from conversation 1
thread_2 = await client.threads.create()
async for chunk in client.runs.stream(
 thread_2["thread_id"],
 "agent",
 input={"messages": [{"role": "human", "content": "What packages are installed?"}]},
 stream_mode="updates",
):
 print(chunk.data)

```
All threads share one sandbox. This is useful when the sandbox has state that’s expensive to recreate, such as a cloned repo, installed dependencies, or build artifacts. Any conversation on the same assistant picks up where the last one left off without repeating setup.client.py
```
from langgraph_sdk import get_client

client = get_client(url="<DEPLOYMENT_URL>", api_key="<LANGSMITH_API_KEY>")

# Conversation 1: clone and set up the project
thread_1 = await client.threads.create()
async for chunk in client.runs.stream(
 thread_1["thread_id"],
 "agent",
 input={"messages": [{"role": "human", "content": "Clone https://github.com/org/repo and install dependencies"}]},
 stream_mode="updates",
):
 print(chunk.data)

# Conversation 2: repo and dependencies are still there
thread_2 = await client.threads.create()
async for chunk in client.runs.stream(
 thread_2["thread_id"],
 "agent",
 input={"messages": [{"role": "human", "content": "Run the test suite and fix any failures"}]},
 stream_mode="updates",
):
 print(chunk.data)

```

#### [​](#file-transfers)File transfers

Sandboxes are isolated containers, so your application code can’t directly access files inside them. Use `upload_files()` and `download_files()` to move data across the sandbox boundary:

{' ' * (self.list_depth - 1)}- **Seed the sandbox before the agent runs**: upload user files, [skill](/oss/python/deepagents/skills) scripts, configuration, or [persistent memories](/oss/python/deepagents/memory) so the agent has what it needs from the start

{' ' * (self.list_depth - 1)}- **Retrieve results after the agent finishes**: download generated artifacts (reports, plots, exports) and sync updated memories back for future conversations

For provider-specific file transfer examples, see [working with files](/oss/python/deepagents/sandboxes#working-with-files). For provider setup, security, and lifecycle patterns, see the full [sandboxes guide](/oss/python/deepagents/sandboxes).

Example: syncing skills and memories with custom middleware[Skill](/oss/python/deepagents/skills) scripts that the agent needs to execute must be uploaded into the sandbox before the agent runs. You may also want to sync [memories](/oss/python/deepagents/memory) so the agent can read and update them inside the container. Use [custom middleware](/oss/python/langchain/middleware/custom) with `before_agent` and `after_agent` hooks to move files across the sandbox boundary:agent.py
```
from deepagents import create_deep_agent
from langchain.agents.middleware import AgentMiddleware, AgentState
from langgraph.runtime import Runtime

def _safe_filename(key: str) -> str:
 """Reject keys that contain path traversal or glob characters."""
 name = key.split("/")[-1]
 if ".." in name or any(c in name for c in ("*", "?")):
 raise ValueError(f"Invalid key: {key}")
 return name

class SandboxSyncMiddleware(AgentMiddleware):
 """Sync skills and memories between the store and the sandbox."""

 def __init__(self, backend: CompositeBackend):
 super().__init__()
 self.backend = backend

 async def abefore_agent(self, state: AgentState, runtime: Runtime) -> None:
 """Upload skill scripts and memories into the sandbox."""
 user_id = runtime.server_info.user.identity 
 store = runtime.store
 files = []
 for item in await store.asearch(("skills", user_id)):
 name = _safe_filename(item.key)
 files.append((f"/skills/{name}", item.value["content"].encode()))
 for item in await store.asearch(("memories", user_id)):
 name = _safe_filename(item.key)
 files.append((f"/memories/{name}", item.value["content"].encode()))
 if files:
 await self.backend.upload_files(files)

 async def aafter_agent(self, state: AgentState, runtime: Runtime) -> None:
 """Sync updated memories back to the store."""
 user_id = runtime.server_info.user.identity 
 store = runtime.store
 items = await store.asearch(("memories", user_id))
 results = await self.backend.download_files(
 [f"/memories/{item.key}" for item in items]
 )
 for result in results:
 if result.content is not None:
 await store.aput(
 ("memories", user_id),
 result.path.split("/")[-1],
 {"content": result.content.decode()},
 )

backend = CompositeBackend(
 default=DaytonaSandbox(sandbox=sandbox),
 routes={
 "/skills/": StoreBackend(
 rt,
 namespace=lambda rt: ("skills", rt.server_info.user.identity),
 ),
 "/memories/": StoreBackend(
 rt,
 namespace=lambda rt: ("memories", rt.server_info.user.identity),
 ),
 },
)

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=backend,
 middleware=[SandboxSyncMiddleware(backend)],
)

```

#### [​](#managing-secrets)Managing secrets

Sandboxes are isolated containers, so environment variables from your host aren’t available inside them. There are two ways to provide API keys and other secrets to sandbox code:
**Auth proxy (recommended).** The [sandbox auth proxy](/langsmith/sandbox-auth-proxy) intercepts outbound requests from the sandbox and injects authentication headers automatically. Sandbox code calls external APIs normally, and the proxy adds the correct credentials based on the destination host. This means API keys never appear in sandbox code, environment variables, or logs.

![The sandbox auth proxy injects credentials into outbound requests so secrets never enter the sandbox](https://mintcdn.com/langchain-5e9cc07a/9oyV6nbtSbBRfaE1/oss/images/deepagents/production/sandbox-auth-proxy.png?fit=max&auto=format&n=9oyV6nbtSbBRfaE1&q=85&s=632c4a493f1d5928e41c6865ab86d1da)

```
{
 "proxy_config": {
 "rules": [
 {
 "name": "openai-api",
 "match_hosts": ["api.openai.com"],
 "inject_headers": {
 "Authorization": "Bearer ${OPENAI_API_KEY}"
 }
 },
 {
 "name": "anthropic-api",
 "match_hosts": ["api.anthropic.com"],
 "inject_headers": {
 "x-api-key": "${ANTHROPIC_API_KEY}"
 }
 }
 ]
 }
}

```

The `${SECRET_KEY}` references resolve against secrets stored in your LangSmith [workspace settings](/langsmith/set-up-hierarchy#configure-workspace-settings). Configure secrets there before creating a template that references them.
**Workspace secrets.** For API keys that don’t need proxy-based injection (e.g., keys used by the agent server itself, not sandbox code), store them as [workspace secrets](/langsmith/set-up-hierarchy#configure-workspace-settings) in LangSmith. These are available as environment variables at runtime for all agents in the workspace.
Avoid passing secrets into sandboxes via environment variables or file uploads. Agents can read any accessible file or environment variable inside the sandbox, including credentials. The auth proxy keeps secrets out of the sandbox entirely.

## [​](#guardrails)Guardrails

Agents in production run autonomously, which means they can loop indefinitely, hit rate limits, or process user data that contains sensitive information. Deep Agents provide two layers of protection:

{' ' * (self.list_depth - 1)}- **[Permissions](/oss/python/deepagents/permissions)**: declarative allow/deny rules that control which files and directories the agent can read or write. Use permissions to isolate the agent to a working directory, protect sensitive files, or enforce read-only memory.

{' ' * (self.list_depth - 1)}- **[Middleware](/oss/python/langchain/middleware/built-in)**: hooks that wrap model and tool calls for rate limiting, error handling, and data privacy.

![Middleware hooks—before_model, wrap_model_call, wrap_tool_call, after_model—wrap the agent loop so policies run deterministically around every relevant step](https://mintcdn.com/langchain-5e9cc07a/9oyV6nbtSbBRfaE1/oss/images/deepagents/production/middleware-lifecycle.png?fit=max&auto=format&n=9oyV6nbtSbBRfaE1&q=85&s=0d30b34aba2b829a1b763b975cfb2817)

### [​](#rate-limiting)Rate limiting

Rate limiting here refers to capping the agent’s own LLM and tool usage within a run, not API gateway rate limiting for incoming requests.
Without limits, a confused agent can burn through your LLM API budget in minutes by looping on the same tool call or making hundreds of model calls. Set caps on both model calls and tool executions per run:

```
from deepagents import create_deep_agent
from langchain.agents.middleware import ModelCallLimitMiddleware, ToolCallLimitMiddleware

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 middleware=[
 ModelCallLimitMiddleware(run_limit=50),
 ToolCallLimitMiddleware(run_limit=200),
 ],
)

```

Use `run_limit` to cap calls within a single invocation (resets each turn). Use `thread_limit` to cap calls across an entire conversation (requires a checkpointer). See [ModelCallLimitMiddleware](https://reference.langchain.com/python/langchain/agents/middleware/model_call_limit/ModelCallLimitMiddleware) and [ToolCallLimitMiddleware](https://reference.langchain.com/python/langchain/agents/middleware/tool_call_limit/ToolCallLimitMiddleware) for the full configuration.

### [​](#handling-errors)Handling errors

Not all errors should be handled the same way. Transient failures (network timeouts, rate limits) should be retried automatically. Errors the LLM can recover from (bad tool output, parsing failures) should be fed back to the model. Errors that need human input should pause the agent. For the full breakdown with code examples, see [Handle errors appropriately](/oss/python/langgraph/thinking-in-langgraph#handle-errors-appropriately).
Middleware handles the transient case. Model calls and tool calls each have their own retry middleware with exponential backoff. If your primary model provider goes down entirely, the fallback middleware switches to an alternative:

```
from deepagents import create_deep_agent
from langchain.agents.middleware import (
 ModelFallbackMiddleware,
 ModelRetryMiddleware,
 ToolRetryMiddleware,
)

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 middleware=[
 # Retry model calls on rate limits, timeouts, and 5xx errors
 ModelRetryMiddleware(max_retries=3, backoff_factor=2.0, initial_delay=1.0),
 # If the primary model is fully down, fall back to an alternative
 ModelFallbackMiddleware("gpt-5.4"),
 # Retry specific tools that hit external APIs (not all tools)
 ToolRetryMiddleware(
 max_retries=2,
 tools=["search", "fetch_url"],
 retry_on=(TimeoutError, ConnectionError),
 ),
 ],
)

```

Scope [ToolRetryMiddleware](https://reference.langchain.com/python/langchain/agents/middleware/tool_retry/ToolRetryMiddleware) to specific tools rather than retrying everything. A filesystem `read_file` that fails won’t benefit from a retry, but a web search that times out probably will. See [ModelRetryMiddleware](https://reference.langchain.com/python/langchain/agents/middleware/model_retry/ModelRetryMiddleware) and [ModelFallbackMiddleware](https://reference.langchain.com/python/langchain/agents/middleware/model_fallback/ModelFallbackMiddleware) for the full configuration.

### [​](#data-privacy)Data privacy

If your agent processes user input that might contain emails, credit card numbers, or other PII, you can detect and handle it before it reaches the model or gets stored in logs:

```
from deepagents import create_deep_agent
from langchain.agents.middleware import PIIMiddleware

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 middleware=[
 PIIMiddleware("email", strategy="redact", apply_to_input=True),
 PIIMiddleware("credit_card", strategy="mask", apply_to_input=True),
 ],
)

```

Strategies include `redact` (replace with `[REDACTED_EMAIL]`), `mask` (partial masking like `****-****-****-1234`), `hash` (deterministic hash), and `block` (raise an error). You can also write custom detectors for domain-specific patterns.
See [PIIMiddleware](https://reference.langchain.com/python/langchain/agents/middleware/pii/PIIMiddleware) for the full configuration.
For the complete list of available middleware, see [prebuilt middleware](/oss/python/langchain/middleware/built-in).

## [​](#frontend)Frontend

Deep Agents use [`useStream`](/oss/python/langchain/frontend/overview) to connect your UI to the agent backend. `useStream` is a frontend hook (available for React, Vue, Svelte, and Angular) that streams messages, subagent progress, and custom state from your agent in real time.
Locally, `useStream` points at `http://localhost:2024`. In production, point it at your [LangSmith Deployment](/langsmith/deployment) and configure reconnection so users don’t lose progress if their connection drops.

```
import { useStream } from "@langchain/react";

function App() {
 const stream = useStream<typeof agent>({
 apiUrl: "https://your-deployment.langsmith.dev",
 assistantId: "agent",
 reconnectOnMount: true, // Resume stream after page refresh or navigation
 fetchStateHistory: true, // Load full thread history on mount
 });
}

```

`reconnectOnMount` picks up an in-progress run automatically. If a user refreshes while the agent is working, they’ll see it continue rather than a blank screen. `fetchStateHistory` loads the full conversation history for the thread, so returning users see previous messages.
For deep agent workflows that spawn many subagents, set a high `recursionLimit` when submitting to avoid cutting off long-running executions:

```
stream.submit(
 { messages: [{ type: "human", content: text }] },
 {
 streamSubgraphs: true,
 config: { recursionLimit: 10000 },
 },
);

```

For UI patterns specific to deep agents, such as subagent cards, todo lists, and custom state rendering, see the [frontend guide](/oss/python/deepagents/frontend/overview).

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/going-to-production.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Managed Deep AgentsPrevious](/langsmith/deploy-managed-deep-agent)[Harness capabilitiesNext](/oss/python/deepagents/harness)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
