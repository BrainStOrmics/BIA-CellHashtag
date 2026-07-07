# Runtime

[Advanced usage](/oss/python/langchain/guardrails)

# Runtime
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.

## [​](#overview)Overview

LangChain’s [`create_agent`](https://reference.langchain.com/python/langchain/agents/factory/create_agent) runs on LangGraph’s runtime under the hood.
LangGraph exposes a [`Runtime`](https://reference.langchain.com/python/langgraph/runtime/Runtime) object with the following information:

{' ' * (self.list_depth - 1)}- **Context**: static information like user id, db connections, or other dependencies for an agent invocation

{' ' * (self.list_depth - 1)}- **Store**: a [BaseStore](https://reference.langchain.com/python/langchain-core/stores/BaseStore) instance used for [long-term memory](/oss/python/langchain/long-term-memory)

{' ' * (self.list_depth - 1)}- **Stream writer**: an object used for streaming information via the `"custom"` stream mode

{' ' * (self.list_depth - 1)}- **Execution info**: identity and retry information for the current execution (thread ID, run ID, attempt number)

{' ' * (self.list_depth - 1)}- **Server info**: server-specific metadata when running on LangGraph Server (assistant ID, graph ID, authenticated user)

Runtime context provides **dependency injection** for your tools and middleware. Instead of hardcoding values or using global state, you can inject runtime dependencies (like database connections, user IDs, or configuration) when invoking your agent. This makes your tools more testable, reusable, and flexible.
You can access the runtime information within [tools](#inside-tools) and [middleware](#inside-middleware).

## [​](#access)Access

When creating an agent with [`create_agent`](https://reference.langchain.com/python/langchain/agents/factory/create_agent), you can specify a `context_schema` to define the structure of the `context` stored in the agent [`Runtime`](https://reference.langchain.com/python/langgraph/runtime/Runtime).
When invoking the agent, pass the `context` argument with the relevant configuration for the run:

```
from dataclasses import dataclass

from langchain.agents import create_agent

@dataclass
class Context:
 user_name: str

agent = create_agent(
 model="gpt-5-nano",
 tools=[...],
 context_schema=Context 
)

agent.invoke(
 {"messages": [{"role": "user", "content": "What's my name?"}]},
 context=Context(user_name="John Smith")
)

```

### [​](#inside-tools)Inside tools

You can access the runtime information inside tools to:

{' ' * (self.list_depth - 1)}- Access the context

{' ' * (self.list_depth - 1)}- Read or write long-term memory

{' ' * (self.list_depth - 1)}- Write to the [custom stream](/oss/python/langchain/streaming#custom-updates) (ex, tool progress / updates)

Use the `ToolRuntime` parameter to access the [`Runtime`](https://reference.langchain.com/python/langgraph/runtime/Runtime) object inside a tool.

```
from dataclasses import dataclass
from langchain.tools import tool, ToolRuntime 

@dataclass
class Context:
 user_id: str

@tool
def fetch_user_email_preferences(runtime: ToolRuntime[Context]) -> str:
 """Fetch the user's email preferences from the store."""
 user_id = runtime.context.user_id 

 preferences: str = "The user prefers you to write a brief and polite email."
 if runtime.store:
 if memory := runtime.store.get(("users",), user_id):
 preferences = memory.value["preferences"]

 return preferences

```

### [​](#execution-info-and-server-info-inside-tools)Execution info and server info inside tools

Access execution identity (thread ID, run ID) via `runtime.execution_info`, and server-specific metadata (assistant ID, authenticated user) via `runtime.server_info` when running on LangGraph Server:

```
from langchain.tools import tool, ToolRuntime

@tool
def context_aware_tool(runtime: ToolRuntime) -> str:
 """A tool that uses execution and server info."""
 # Access thread and run IDs
 info = runtime.execution_info
 print(f"Thread: {info.thread_id}, Run: {info.run_id}")

 # Access server info (only available on LangGraph Server)
 server = runtime.server_info
 if server is not None:
 print(f"Assistant: {server.assistant_id}")
 if server.user is not None:
 print(f"User: {server.user.identity}")

 return "done"

```

`server_info` is `None` when not running on LangGraph Server (e.g., during local development).
Requires `deepagents>=0.5.0` (or `langgraph>=1.1.5`) for `runtime.execution_info` and `runtime.server_info`.

### [​](#inside-middleware)Inside middleware

You can access runtime information in middleware to create dynamic prompts, modify messages, or control agent behavior based on user context.
Use the `Runtime` parameter to access the [`Runtime`](https://reference.langchain.com/python/langgraph/runtime/Runtime) object inside [node-style hooks](/oss/python/langchain/middleware/custom#node-style-hooks). For [wrap-style hooks](/oss/python/langchain/middleware/custom#wrap-style-hooks), the `Runtime` object is available inside the [`ModelRequest`](https://reference.langchain.com/python/langchain/agents/middleware/types/ModelRequest) parameter.

```
from dataclasses import dataclass

from langchain.messages import AnyMessage
from langchain.agents import create_agent, AgentState
from langchain.agents.middleware import dynamic_prompt, ModelRequest, before_model, after_model
from langgraph.runtime import Runtime

@dataclass
class Context:
 user_name: str

# Dynamic prompts
@dynamic_prompt
def dynamic_system_prompt(request: ModelRequest) -> str:
 user_name = request.runtime.context.user_name 
 system_prompt = f"You are a helpful assistant. Address the user as {user_name}."
 return system_prompt

# Before model hook
@before_model
def log_before_model(state: AgentState, runtime: Runtime[Context]) -> dict | None:
 print(f"Processing request for user: {runtime.context.user_name}")
 return None

# After model hook
@after_model
def log_after_model(state: AgentState, runtime: Runtime[Context]) -> dict | None:
 print(f"Completed request for user: {runtime.context.user_name}")
 return None

agent = create_agent(
 model="gpt-5-nano",
 tools=[...],
 middleware=[dynamic_system_prompt, log_before_model, log_after_model],
 context_schema=Context
)

agent.invoke(
 {"messages": [{"role": "user", "content": "What's my name?"}]},
 context=Context(user_name="John Smith")
)

```

### [​](#execution-info-and-server-info-inside-middleware)Execution info and server info inside middleware

Middleware hooks can also access `runtime.execution_info` and `runtime.server_info`:

```
from langchain.agents import AgentState
from langchain.agents.middleware import before_model
from langgraph.runtime import Runtime

@before_model
def auth_gate(state: AgentState, runtime: Runtime) -> dict | None:
 """Block unauthenticated users when running on LangGraph Server."""
 server = runtime.server_info
 if server is not None and server.user is None:
 raise ValueError("Authentication required")
 print(f"Thread: {runtime.execution_info.thread_id}")
 return None

```

Requires `deepagents>=0.5.0` (or `langgraph>=1.1.5`).

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/runtime.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[GuardrailsPrevious](/oss/python/langchain/guardrails)[Context engineering in agentsNext](/oss/python/langchain/context-engineering)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
