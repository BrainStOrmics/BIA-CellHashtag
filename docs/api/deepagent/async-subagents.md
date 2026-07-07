# Async subagents

[Core capabilities](/oss/python/deepagents/harness)

# Async subagents
Copy page

Launch background subagents that run concurrently while the supervisor continues interacting with the userCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Async subagents let a supervisor agent launch background tasks that return immediately, so the supervisor can continue interacting with the user while subagents work concurrently. The supervisor can check progress, send follow-up instructions, or cancel tasks at any point.
This builds on [subagents](/oss/python/deepagents/subagents), which run synchronously and block the supervisor until completion. Use async subagents when tasks are long-running, parallelizable, or need mid-flight steering.
Async subagents are a preview feature available in `deepagents` 0.5.0. Preview features are under active development and APIs may change.

Async subagents communicate with any server that implements the [Agent Protocol](https://github.com/langchain-ai/agent-protocol). You can use [LangSmith Deployments](/langsmith/deployment), or self-host any Agent Protocol-compatible server. Each subagent runs independently of the supervisor, which controls them through the SDK to launch, check, update, and cancel.

## [​](#when-to-use-async-subagents)When to use async subagents

| 
 | Dimension | Sync subagents | Async subagents
 | **Execution model** | Supervisor blocks until subagent completes | Returns job ID immediately; supervisor continues
 | **Concurrency** | Parallel but blocking | Parallel and non-blocking
 | **Mid-task updates** | Not possible | Send follow-up instructions via `update_async_task`
 | **Cancellation** | Not possible | Cancel running tasks via `cancel_async_task`
 | **Statefulness** | Stateless — no persistent state between invocations | Stateful — maintains state on its own thread across interactions
 | **Best for** | Tasks where the agent should wait for results before continuing | Long-running, complex tasks managed interactively in a chat

## [​](#configure-async-subagents)Configure async subagents

Define async subagents as a list of [`AsyncSubAgent`](https://reference.langchain.com/python/deepagents/middleware/async_subagents/AsyncSubAgent) specs, each pointing to an Agent Protocol server:

```
from deepagents import AsyncSubAgent, create_deep_agent

async_subagents = [
 AsyncSubAgent(
 name="researcher",
 description="Research agent for information gathering and synthesis",
 graph_id="researcher",
 # No url → ASGI transport (co-deployed in the same deployment)
 ),
 AsyncSubAgent(
 name="coder",
 description="Coding agent for code generation and review",
 graph_id="coder",
 # url="https://coder-deployment.langsmith.dev" # Optional: HTTP transport for remote
 ),
]

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 subagents=async_subagents,
)

```

| 
 | Field | Type | Description
 | `name` | `str` | Required. Unique identifier. The supervisor uses this when launching tasks.
 | `description` | `str` | Required. What this subagent does. The supervisor uses this to decide which agent to delegate to.
 | `graph_id` | `str` | Required. The graph ID (or assistant ID) on the Agent Protocol server. For LangGraph-based deployments, this must match a graph registered in `langgraph.json`.
 | `url` | `str` | Optional. When omitted, uses ASGI transport (in-process). When set, uses HTTP transport to a remote Agent Protocol server.
 | `headers` | `dict[str, str]` | Optional. Additional headers for requests to the remote server. Use for custom authentication with self-hosted Agent Protocol servers.
For LangGraph-based deployments, register all graphs in the same `langgraph.json` for co-deployed setups:

```
{
 "graphs": {
 "supervisor": "./src/supervisor.py:graph",
 "researcher": "./src/researcher.py:graph",
 "coder": "./src/coder.py:graph"
 }
}

```

## [​](#use-the-async-subagent-tools)Use the async subagent tools

The [`AsyncSubAgentMiddleware`](https://reference.langchain.com/python/deepagents/middleware/async_subagents/AsyncSubAgentMiddleware) gives the supervisor five tools:

| 
 | Tool | Purpose | Returns
 | `start_async_task` | Start a new background task | Task ID (immediately)
 | `check_async_task` | Get current status and result of a task | Status + result (if complete)
 | `update_async_task` | Send new instructions to a running task | Confirmation + updated status
 | `cancel_async_task` | Stop a running task | Confirmation
 | `list_async_tasks` | List all tracked tasks with live statuses | Summary of all tasks
The supervisor’s LLM calls these tools like any other tool. The middleware handles thread creation, run management, and state persistence automatically.

### [​](#understand-the-lifecycle)Understand the lifecycle

A typical interaction follows this sequence:

{' ' * (self.list_depth - 1)}- **Launch** creates a new thread on the server, starts a run with the task description as input, and returns the thread ID as the task ID. The supervisor reports this ID to the user and does not poll for completion.

{' ' * (self.list_depth - 1)}- **Check** fetches the current run status. If the run succeeded, it retrieves the thread state to extract the subagent’s final output. If still running, it reports that to the user.

{' ' * (self.list_depth - 1)}- **Update** creates a new run on the same thread with an interrupt multitask strategy. The previous run is interrupted, and the subagent restarts with the full conversation history plus the new instructions. The task ID stays the same.

{' ' * (self.list_depth - 1)}- **Cancel** calls `runs.cancel()` on the server and marks the task as `"cancelled"`.

{' ' * (self.list_depth - 1)}- **List** iterates over all tracked tasks. For non-terminal tasks, it fetches live status from the server in parallel. Terminal statuses (`success`, `error`, `cancelled`) are returned from cache.

## [​](#understand-state-management)Understand state management

Task metadata is stored in a dedicated state channel (`async_tasks`) on the supervisor’s graph, separate from the message history. This is critical because deep agents [compact their message history](/oss/python/deepagents/context-engineering#summarization) when the context window fills up. If task IDs were only in tool messages, they would be lost during compaction. The dedicated channel ensures the supervisor can always recall its tasks through `list_async_tasks`, even after multiple rounds of summarization.
Each tracked task records the task ID, agent name, thread ID, run ID, status, and timestamps (`created_at`, `last_checked_at`, `last_updated_at`).

## [​](#choose-a-transport)Choose a transport

### [​](#asgi-transport-co-deployed)ASGI transport (co-deployed)

When a subagent spec omits the `url` field, the LangGraph SDK uses ASGI transport — SDK calls are routed through in-process function calls rather than HTTP. For LangGraph-based deployments, this requires both graphs to be registered in the same `langgraph.json`.
ASGI transport eliminates network latency and requires no additional auth configuration. The subagent still runs as a separate thread with its own state. This is the recommended default.

### [​](#http-transport-remote)HTTP transport (remote)

Add a `url` field to switch to HTTP transport, where SDK calls go over the network to a remote Agent Protocol server:

```
AsyncSubAgent(
 name="researcher",
 description="Research agent",
 graph_id="researcher",
 url="https://my-research-deployment.langsmith.dev",
)

```

For LangGraph deployments, authentication is handled by the LangGraph SDK using `LANGSMITH_API_KEY` (or `LANGGRAPH_API_KEY`) from environment variables. Self-hosted Agent Protocol servers may use a different authentication mechanism.
Use HTTP transport when subagents need independent scaling, different resource profiles, or are maintained by a different team.

## [​](#choose-a-deployment-topology)Choose a deployment topology

### [​](#single-deployment)Single deployment

A single deployment means all agents are co-deployed on the same server using ASGI transport. For LangGraph-based deployments, register all graphs in one `langgraph.json`. This is the recommended starting point — one server to manage, zero network latency between agents.

### [​](#split-deployment)Split deployment

Supervisor on one server, subagents on another via HTTP transport. Use when subagents need different compute profiles or independent scaling.

### [​](#hybrid)Hybrid

In a split deployment, you have the supervisor on one server and subagents on another via HTTP transport. Use when subagents need different compute profiles or independent scaling.

### [​](#hybrid-2)Hybrid

In a hybrid deployment, some subagents are co-deployed via ASGI, others remote via HTTP:

```
async_subagents = [
 AsyncSubAgent(
 name="researcher",
 description="Research agent",
 graph_id="researcher",
 # No url → ASGI (co-deployed)
 ),
 AsyncSubAgent(
 name="coder",
 description="Coding agent",
 graph_id="coder",
 url="https://coder-deployment.langsmith.dev",
 # url present → HTTP (remote)
 ),
]

```

## [​](#best-practices)Best practices

### [​](#size-the-worker-pool-for-local-development)Size the worker pool for local development

When running locally with `langgraph dev`, increase the worker pool to accommodate concurrent subagent runs. Each active run occupies a worker slot. A supervisor with 3 concurrent subagent tasks requires 4 slots (1 supervisor + 3 subagents). Under-provisioning causes launches to queue.

```
langgraph dev --n-jobs-per-worker 10

```

### [​](#write-clear-subagent-descriptions)Write clear subagent descriptions

The supervisor uses descriptions to decide which subagent to launch. Be specific and action-oriented:

```
# Good
AsyncSubAgent(
 name="researcher",
 description="Conducts in-depth research using web search. Use for questions requiring multiple searches and synthesis.",
 graph_id="researcher",
)

# Bad
AsyncSubAgent(
 name="helper",
 description="helps with stuff",
 graph_id="helper",
)

```

### [​](#trace-with-thread-ids)Trace with thread IDs

When using LangGraph-based deployments, every async subagent run is a standard LangGraph run, fully visible in LangSmith. The supervisor’s trace shows tool calls for `launch`, `check`, `update`, `cancel`, and `list`. Each subagent run appears as a separate trace, linked by thread ID. Use the thread ID (task ID) to correlate supervisor orchestration traces with subagent execution traces.

## [​](#troubleshooting)Troubleshooting

### [​](#supervisor-polls-immediately-after-launch)Supervisor polls immediately after launch

**Problem**: The supervisor calls `check` in a loop right after launching, turning async execution into blocking.
**Solution**: The middleware injects system prompt rules to prevent this. If polling persists, reinforce the behavior in your supervisor’s system prompt:

```
agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 system_prompt="""...your instructions...

 After launching an async subagent, ALWAYS return control to the user.
 Never call check_async_task immediately after launch.""",
 subagents=async_subagents,
)

```

### [​](#supervisor-reports-stale-status)Supervisor reports stale status

**Problem**: The supervisor references a task status from earlier in conversation history instead of making a fresh `check` call.
**Solution**: The middleware prompt instructs the model that “task statuses in conversation history are always stale.” If this still occurs, add explicit instructions to always call `check` or `list` before reporting status.

### [​](#task-id-lookup-failures)Task ID lookup failures

**Problem**: The supervisor truncates or reformats the task ID, causing `check` or `cancel` to fail.
**Solution**: The middleware prompt instructs the model to always use the full task ID. If truncation persists, this is typically a model-specific issue — try a different model or add “always show the full task_id, never truncate or abbreviate it” to your system prompt.

### [​](#subagent-launches-queue-instead-of-running)Subagent launches queue instead of running

**Problem**: Launching a subagent hangs or takes a long time to start.
**Solution**: The worker pool is likely exhausted. Increase the pool size with `--n-jobs-per-worker`. See [Size the worker pool](#size-the-worker-pool-for-local-development).

## [​](#reference-implementation)Reference implementation

The [async-deep-agents](https://github.com/langchain-ai/async-deep-agents) repository contains working examples in both Python and TypeScript that deploy to LangSmith Deployments. It demonstrates a supervisor with researcher and coder subagents running as background tasks.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/async-subagents.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[SubagentsPrevious](/oss/python/deepagents/subagents)[Human-in-the-loopNext](/oss/python/deepagents/human-in-the-loop)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
