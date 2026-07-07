# Context engineering in Deep Agents

[Core capabilities](/oss/python/deepagents/harness)

# Context engineering in Deep Agents
Copy page

Control what context your deep agent has access to and how it is managed across long-running tasksCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Context engineering is providing the right information and tools in the right format so your deep agent can accomplish tasks reliably.
Deep agents have access to several kinds of context.
Some sources are provided to the agent at startup; others become available during runtime, such as user input.
Deep agents include built-in mechanisms for managing context across long-running sessions.
This page provides an overview of the different kinds of context your deep agent has access to and manages.
New to context engineering? See the [conceptual overview](/oss/python/concepts/context) for the different types of context and when to use them.

## [​](#types-of-context)Types of context

| 
 | Context Type | What You Control | Scope
 | **[Input context](#input-context)** | What goes into the agent’s prompt at startup (system prompt, memory, skills) | Static, applied each run
 | **[Runtime context](#runtime-context)** | Static configuration passed at invoke time (user metadata, API keys, connections) | Per run, propagates to subagents
 | **[Context compression](#context-compression)** | Built-in offloading and summarization to keep context within window limits | Automatic, when limits approached
 | **[Context isolation](#context-isolation-with-subagents)** | Use subagents to quarantine heavy work, returning only results to the main agent | Per subagent, when delegated
 | **[Long-term memory](#long-term-memory)** | Persistent storage across threads using the virtual filesystem | Persistent across conversations

## [​](#input-context)Input context

Input context is information provided to your deep agent at startup that becomes part of its system prompt. The final prompt consists of several sources:

## System prompt
Custom instructions you provide plus built-in agent guidance.

## Memory
Persistent `AGENTS.md` files always loaded when configured.

## Skills
On-demand capabilities loaded when relevant (progressive disclosure).

## Tool prompts
Instructions for using built-in tools or custom tools.

### [​](#system-prompt)System prompt

Your custom system prompt is prepended to the built-in system prompt, which includes guidance for planning, filesystem tools, and subagents. Use it to define the agent’s role, behavior, and knowledge:

```
from deepagents import create_deep_agent

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 system_prompt=(
 "You are a research assistant specializing in scientific literature. "
 "Always cite sources. Use subagents for parallel research on different topics."
 ),
)

```

The `system_prompt` parameter is static which means it does not change per invocation.
For some use cases you may want a dynamic prompt: for example, to tell the model “You have admin access” vs “You have read-only access,” or to inject user preferences like “User prefers concise responses” from [long-term memory](#long-term-memory).
If your prompt depends on context or `runtime.store`, use `@dynamic_prompt` to build context-aware instructions. Your middleware can read `request.runtime.context` and `request.runtime.store`.
See [Customization](/oss/python/deepagents/customization#middleware) for adding [custom middleware](/oss/python/langchain/middleware) and the [LangChain context engineering](/oss/python/langchain/context-engineering#system-prompt) guide for examples.
You do **not** need middleware when tools alone use context or `runtime.store`; tools receive the [ToolRuntime](https://reference.langchain.com/python/langchain/tools/#langchain.tools.ToolRuntime) object (including `runtime.context` and `runtime.store`) directly. Add middleware only when tools should be packaged with an update to the system prompt.
To adjust the assembled system prompt for a specific provider or model, use a [harness profile](/oss/python/deepagents/profiles#harness-profiles): `base_system_prompt` replaces the base prompt outright, and `system_prompt_suffix` appends to it.

### [​](#memory)Memory

Memory files ([`AGENTS.md`](https://agents.md/)) provide persistent context that is **always loaded** into the system prompt. Use memory for project conventions, user preferences, and critical guidelines that should apply to every conversation:

```
agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 memory=["/project/AGENTS.md", "~/.deepagents/preferences.md"],
)

```

Unlike skills, memory is always injected—there is no progressive disclosure. Keep memory minimal to avoid context overload; use [skills](/oss/python/deepagents/skills) for detailed workflows and domain-specific content. See [Memory](/oss/python/deepagents/customization#memory) for configuration details.

### [​](#skills)Skills

Skills provide **on-demand** capabilities. The agent reads frontmatter from each `SKILL.md` at startup, then loads full skill content only when it determines the skill is relevant. This reduces token usage while still providing specialized workflows:

```
agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 skills=["/skills/research/", "/skills/web-search/"],
)

```

Keep each skill focused on a single workflow or domain; broad or overlapping skills dilute relevance and bloat context when loaded. Within a skill, keep the main content concise and move detailed reference material to separate files that are referenced in the skill file. Put always-relevant conventions in [memory](#memory). See [Skills](/oss/python/deepagents/skills) for authoring and configuration.

### [​](#tool-prompts)Tool prompts

[Tool](/oss/python/langchain/tools) prompts are instructions that shape how the model uses tools. All tools expose metadata the model sees in its prompt—typically a schema and a description. Tools you pass via the `tools` parameter surface that tool metadata (schema and descriptions) to the model. A deep agent’s built-in tools are packaged in middleware and typically also update the system prompt with more guidance for those tools.
**Built-in tools** – Middleware that adds harness capabilities (planning, filesystem, subagents) automatically appends tool-specific instructions to the system prompt, creating tool prompts that explain how to use those tools effectively:

{' ' * (self.list_depth - 1)}- Planning prompt – Instructions for `write_todos` to maintain a structured task list

{' ' * (self.list_depth - 1)}- Filesystem prompt – Documentation for `ls`, `read_file`, `write_file`, `edit_file`, `glob`, `grep` (and `execute` when using a sandbox backend)

{' ' * (self.list_depth - 1)}- Subagent prompt – Guidance for delegating work with the `task` tool

{' ' * (self.list_depth - 1)}- Human-in-the-loop prompt – Usage for pausing at specified tool calls (when `interrupt_on` is set)

{' ' * (self.list_depth - 1)}- Local context prompt – Current directory and project info (CLI only)

**Tools you provide** – Tools passed via the `tools` parameter get their descriptions (from the tool schema) sent to the model. You can also add [custom middleware](/oss/python/langchain/middleware) that adds tools and appends its own system prompt instructions.
For tools you provide, make sure to provide a clear name, description, and argument descriptions. These guide the model’s reasoning about when and how to use the tool. Include *when* to use the tool in the description and describe what each argument does.

```
@tool(parse_docstring=True)
def search_orders(
 user_id: str,
 status: str,
 limit: int = 10
) -> str:
 """Search for user orders by status.

 Use this when the user asks about order history or wants to check
 order status. Always filter by the provided status.

 Args:
 user_id: Unique identifier for the user
 status: Order status: 'pending', 'shipped', or 'delivered'
 limit: Maximum number of results to return
 """
 # Implementation here
 ...

```

To override a built-in or user-supplied tool’s description for a specific provider or model, use a [harness profile](/oss/python/deepagents/profiles#harness-profiles)’s `tool_description_overrides` keyed by tool name. `excluded_tools` removes a tool from the visible tool set entirely.
See [Harness](/oss/python/deepagents/harness) for built-in capabilities and [Customization](/oss/python/deepagents/customization#tools) for passing tools directly.

### [​](#complete-system-prompt)Complete system prompt

The deep agent’s system message—the assembled system prompt the model receives at the start of a run—consists of the following parts:

{' ' * (self.list_depth - 1)}- Custom `system_prompt` (if provided)

{' ' * (self.list_depth - 1)}- [Base agent prompt](https://github.com/langchain-ai/deepagents/blob/e18e9dcd0e6edc72c0a4a5b76ae752c4bc539752/libs/deepagents/deepagents/graph.py#L37)

{' ' * (self.list_depth - 1)}- To-do list prompt: Instructions for how to plan with to do lists

{' ' * (self.list_depth - 1)}- Memory prompt: `AGENTS.md` + memory usage guidelines (only when `memory` provided)

{' ' * (self.list_depth - 1)}- Skills prompt: Skills locations + list of skills with frontmatter information + usage (only when skills provided)

{' ' * (self.list_depth - 1)}- Virtual filesystem prompt (filesystem + execute tool docs if applicable)

{' ' * (self.list_depth - 1)}- Subagent prompt: Task tool usage

{' ' * (self.list_depth - 1)}- User-provided middleware prompts (if custom middleware is provided)

{' ' * (self.list_depth - 1)}- Human-in-the-loop prompt (when `interrupt_on` is set)

## [​](#runtime-context)Runtime context

Runtime context is per-run configuration you pass when you invoke the agent. It is not automatically included in the model prompt; the model only sees it if a tool, middleware, or other logic reads it and adds it to messages or the system prompt. Use runtime context for user metadata (IDs, preferences, roles), API keys, database connections, feature flags, or other values your tools and harness need.
Define the shape of that data with `context_schema`: use a `dataclasses.dataclass` or `typing.TypedDict` class. Pass values with the **`context`** argument to `invoke` / `ainvoke`. See [Runtime](/oss/python/langchain/runtime) and [LangGraph runtime context](/oss/python/langgraph/graph-api#runtime-context) for full detail.
Inside tools, read context from the injected [ToolRuntime](https://reference.langchain.com/python/langchain/tools/#langchain.tools.ToolRuntime):

```
from dataclasses import dataclass

from deepagents import create_deep_agent
from langchain.tools import tool, ToolRuntime

@dataclass
class Context:
 user_id: str
 api_key: str

@tool
def fetch_user_data(query: str, runtime: ToolRuntime[Context]) -> str:
 """Fetch data for the current user."""
 user_id = runtime.context.user_id
 return f"Data for user {user_id}: {query}"

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[fetch_user_data],
 context_schema=Context,
)

result = agent.invoke(
 {"messages": [{"role": "user", "content": "Get my recent activity"}]},
 context=Context(user_id="user-123", api_key="sk-..."),
)

```

Runtime context **propagates to all subagents**. When a subagent runs, it receives the same runtime context as the parent. See [Subagents](/oss/python/deepagents/subagents#context-management) for per-subagent context (namespaced keys).

## [​](#context-compression)Context compression

Long-running tasks produce large tool outputs and long conversation history.
Context compression reduces the size of information in an agent’s working memory while preserving details relevant to the task.
The following techniques are the built-in mechanisms to ensure the context passed to LLMs stays within its context window limit:

## Offloading
Large tool inputs and results are stored in the filesystem and replaced with references.

## Summarization
Old messages are compressed into an LLM-generated summary when limits are approached.

### [​](#offloading)Offloading

Deep Agents use the [built-in filesystem tools](/oss/python/deepagents/harness#virtual-filesystem-access) to automatically offload content and to search and retrieve that offloaded content as needed.
Content offloading happens when tool call inputs or results exceed a token threshold (default 20,000):

{' ' * (self.list_depth - 1)}- 
**Tool call inputs exceed 20,000 tokens**: File write and edit operations leave behind tool calls containing the complete file content in the agent’s conversation history.
Since this content is already persisted to the filesystem, it’s often redundant.
As the session context crosses 85% of the model’s available window, deep agents truncate older tool calls, replacing them with a pointer to the file on disk and reducing the size of the active context.

![An example of offloading showing a large input which is saved to disk and the truncated version is used for the tool call](https://mintcdn.com/langchain-5e9cc07a/0G7fpRWZQ2tFN1wL/oss/images/deepagents/offloading-inputs.png?fit=max&auto=format&n=0G7fpRWZQ2tFN1wL&q=85&s=fa18372080684d661965ea6f5ed1edd0)

{' ' * (self.list_depth - 1)}- 
**Tool call results exceed 20,000 tokens**: When this occurs, the deep agent offloads the response to the configured backend and substitutes it with a file path reference and a preview of the first 10 lines. Agents can then re-read or search the content as needed.

![An example of offloading showing a large tool response that is replaced with a message about the location of the offloaded results and the first 10 lines of the result](https://mintcdn.com/langchain-5e9cc07a/0G7fpRWZQ2tFN1wL/oss/images/deepagents/offloading-results.png?fit=max&auto=format&n=0G7fpRWZQ2tFN1wL&q=85&s=11f3da2f37cae63b8aa4c440549f1a67)

### [​](#summarization)Summarization

When the context size crosses the model’s context window limit (for example 85% of `max_input_tokens`), and there is no more context eligible for offloading, the deep agent summarizes the message history.
This process has two components:

{' ' * (self.list_depth - 1)}- **In-context summary**: An LLM generates a structured summary of the conversation including session intent, artifacts created, and next steps—which replaces the full conversation history in the agent’s working memory.

{' ' * (self.list_depth - 1)}- **Filesystem preservation**: The complete, original conversation messages are written to the filesystem as a canonical record.

This dual approach ensures the agent maintains awareness of its goals and progress (via the summary) while preserving the ability to recover specific details when needed (via filesystem search).

![An example of summarization showing an agent's conversation history, where several steps get compacted](https://mintcdn.com/langchain-5e9cc07a/0G7fpRWZQ2tFN1wL/oss/images/deepagents/summarization.png?fit=max&auto=format&n=0G7fpRWZQ2tFN1wL&q=85&s=a8fea59d4365dd688e49ce118e706e76)
**Configuration:**

{' ' * (self.list_depth - 1)}- Triggers at 85% of the model’s `max_input_tokens` from its [model profile](/oss/python/langchain/models#model-profiles)

{' ' * (self.list_depth - 1)}- Keeps 10% of tokens as recent context

{' ' * (self.list_depth - 1)}- Falls back to 170,000-token trigger / 6 messages kept if model profile is unavailable

{' ' * (self.list_depth - 1)}- If any model call raises a standard [ContextOverflowError](https://reference.langchain.com/python/langchain-core/exceptions/ContextOverflowError), the deep agent immediately falls back to summarization and retry with summary + recent preserved messages

{' ' * (self.list_depth - 1)}- Older messages are summarized by the model

[Streaming tokens](/oss/python/deepagents/streaming#llm-tokens) from the agent will generally include tokens generated by the summarization step. You can filter out these tokens using their associated metadata:
```
for chunk in agent.stream(
 {"messages": [...]},
 stream_mode="messages",
 version="v2",
):
 token, metadata = chunk["data"]
 if metadata.get("lc_source") == "summarization":
 continue
 else:
 ...

```

##### Summarization Tool

Deep agents include an optional [tool](/oss/python/langchain/tools) for summarization, enabling agents to trigger summarization at opportune times—such as between tasks—instead of at fixed token intervals.
You can enable this tool by appending it to the middleware list:

```
from deepagents import create_deep_agent
from deepagents.backends import StateBackend
from deepagents.middleware.summarization import (
 create_summarization_tool_middleware,
)

backend = StateBackend # if using default backend

model = "google_genai:gemini-3.5-flash"
agent = create_deep_agent(
 model=model,
 middleware=[
 create_summarization_tool_middleware(model, backend),
 ],
)

```

Enabling this feature does not disable the default summarization action at 85% of the model’s context limit.
See the [`SummarizationToolMiddleware`](https://reference.langchain.com/python/deepagents/middleware/summarization/SummarizationToolMiddleware) API reference for details.

## [​](#context-isolation-with-subagents)Context isolation with subagents

Subagents solve the **context bloat problem**. When the main agent uses tools with large outputs (web search, file reads, database queries), the context window fills quickly. Subagents isolate this work—the main agent receives only the final result, not the dozens of tool calls that produced it. You can also configure each subagent separately from the main agent (for example, model, tools, system prompt, and skills).
**How it works:**

{' ' * (self.list_depth - 1)}- Main agent has a `task` tool to delegate work

{' ' * (self.list_depth - 1)}- Subagent runs with its own fresh context

{' ' * (self.list_depth - 1)}- Subagent executes autonomously until completion

{' ' * (self.list_depth - 1)}- Subagent returns a single final report to the main agent

{' ' * (self.list_depth - 1)}- Main agent’s context stays clean

**Best practices:**

{' ' * (self.list_depth - 1)}- 
**Delegate complex tasks**: Use subagents for multi-step work that would clutter the main agent’s context.

{' ' * (self.list_depth - 1)}- 
**Keep subagent responses concise**: Instruct subagents to return summaries, not raw data:

```
research_subagent = {
 "name": "researcher",
 "description": "Conducts research on a topic",
 "system_prompt": """You are a research assistant.
 IMPORTANT: Return only the essential summary (under 500 words).
 Do NOT include raw search results or detailed tool outputs.""",
 "tools": [web_search],
}

```

{' ' * (self.list_depth - 1)}- 
**Use the filesystem for large data**: Subagents can write results to files; the main agent reads what it needs.

See [Subagents](/oss/python/deepagents/subagents) for configuration and [context management](/oss/python/deepagents/subagents#context-management) for runtime context propagation and per-subagent namespacing.

## [​](#long-term-memory)Long-term memory

When using the default filesystem, your deep agent stores its working memory files in agent state, which only persists within a single thread.
Long-term memory enables your deep agent to persist information across different threads and conversations.
Deep agents can use long-term memory for storing user preferences, accumulated knowledge, research progress, or any information that should persist beyond a single session.
To use long-term memory, you must use a `CompositeBackend` that routes specific paths (typically `/memories/`) to a LangGraph Store, which provides durable cross-thread persistence.
The `CompositeBackend` is a hybrid storage system where some files persist indefinitely while others remain scoped to a single thread.

```
from deepagents import create_deep_agent
from deepagents.backends import CompositeBackend, StateBackend, StoreBackend
from langgraph.store.memory import InMemoryStore

def make_backend(runtime):
 return CompositeBackend(
 default=StateBackend(runtime),
 routes={"/memories/": StoreBackend(runtime)},
 )

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 store=InMemoryStore(),
 backend=make_backend,
 system_prompt="""When users tell you their preferences, save them to
 /memories/user_preferences.txt so you remember them in future conversations.""",
)

```

You do not need to pre-populate `/memories/` with files.
You provide the backend config, store, and system prompt instructions that tell the agent *what* to save and *where*.
For example, you may prompt the agent to store preferences in `/memories/preferences.txt`.
The path starts empty and the agent creates files on demand using its filesystem tools (`write_file`, `edit_file`) when users share information worth remembering.
To pre-seed memories, use the [Store API](/langsmith/custom-store) when deploying on LangSmith.
See [Long-term memory](/oss/python/deepagents/memory) for setup and use cases.

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Start with the right input context** – Keep memory minimal for always-relevant conventions; use focused skills for task-specific capabilities.

{' ' * (self.list_depth - 1)}- **Leverage subagents for heavy work** – Delegate multi-step, output-heavy tasks to keep the main agent’s context clean.

{' ' * (self.list_depth - 1)}- **Adjust subagent outputs in configuration** – If you notice when debugging that subagents generate long output, you can add guidance to the subagent’s `system_prompt` to create summaries and synthesized findings.

{' ' * (self.list_depth - 1)}- **Use the filesystem** – Persist large outputs to files (for example subagent writes or [automatic offloading](#offloading)) so the active context stays small; the model can pull in fragments with `read_file` and `grep` when it needs details.

{' ' * (self.list_depth - 1)}- **Document long-term memory structure** – Tell the agent what lives in `/memories/` and how to use it.

{' ' * (self.list_depth - 1)}- **Pass runtime context for tools** – Use `context` for user metadata, API keys, and other static configuration that tools need.

## [​](#related-resources)Related resources

{' ' * (self.list_depth - 1)}- [Harness](/oss/python/deepagents/harness) – Context management overview, offloading, summarization

{' ' * (self.list_depth - 1)}- [Subagents](/oss/python/deepagents/subagents) – Context isolation, runtime context propagation

{' ' * (self.list_depth - 1)}- [Long-term memory](/oss/python/deepagents/memory) – Cross-thread persistence

{' ' * (self.list_depth - 1)}- [Skills](/oss/python/deepagents/skills) – Progressive disclosure and skill authoring

{' ' * (self.list_depth - 1)}- [Backends](/oss/python/deepagents/backends) – Filesystem backends and CompositeBackend

{' ' * (self.list_depth - 1)}- [Context conceptual overview](/oss/python/concepts/context) – Context types and lifecycle

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/context-engineering.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[ModelsPrevious](/oss/python/deepagents/models)[BackendsNext](/oss/python/deepagents/backends)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
