# Tools

[Core components](/oss/python/langchain/agents)

# Tools
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Tools extend what [agents](/oss/python/langchain/agents) can do—letting them fetch real-time data, execute code, query external databases, and take actions in the world.
Under the hood, tools are callable functions with well-defined inputs and outputs that get passed to a [chat model](/oss/python/langchain/models). The model decides when to invoke a tool based on the conversation context, and what input arguments to provide.
For details on how models handle tool calls, see [Tool calling](/oss/python/langchain/models#tool-calling). Trace tool calls and debug errors with [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-tools)—follow the [tracing quickstart](/langsmith/trace-with-langchain) to get set up.

## [​](#create-tools)Create tools

### [​](#basic-tool-definition)Basic tool definition

The simplest way to create a tool is with the [`@tool`](https://reference.langchain.com/python/langchain-core/tools/convert/tool) decorator. By default, the function’s docstring becomes the tool’s description that helps the model understand when to use it:

```
from langchain.tools import tool

@tool
def search_database(query: str, limit: int = 10) -> str:
 """Search the customer database for records matching the query.

 Args:
 query: Search terms to look for
 limit: Maximum number of results to return
 """
 return f"Found {limit} results for '{query}'"

```

Type hints are **required** as they define the tool’s input schema. The docstring should be informative and concise to help the model understand the tool’s purpose.
**Server-side tool use:** Some chat models feature built-in tools (web search, code interpreters) that are executed server-side. See [Server-side tool use](#server-side-tool-use) for details.
Prefer `snake_case` for tool names (e.g., `web_search` instead of `Web Search`). Some model providers have issues with or reject names containing spaces or special characters with errors. Sticking to alphanumeric characters, underscores, and hyphens helps to improve compatibility across providers.

### [​](#customize-tool-properties)Customize tool properties

#### [​](#custom-tool-name)Custom tool name

By default, the tool name comes from the function name. Override it when you need something more descriptive:

```
@tool("web_search") # Custom name
def search(query: str) -> str:
 """Search the web for information."""
 return f"Results for: {query}"

print(search.name) # web_search

```

#### [​](#custom-tool-description)Custom tool description

Override the auto-generated tool description for clearer model guidance:

```
@tool("calculator", description="Performs arithmetic calculations. Use this for any math problems.")
def calc(expression: str) -> str:
 """Evaluate mathematical expressions."""
 return str(eval(expression))

```

### [​](#advanced-schema-definition)Advanced schema definition

Define complex inputs with Pydantic models or JSON schemas:
Pydantic modelJSON Schema
```
from pydantic import BaseModel, Field
from typing import Literal

class WeatherInput(BaseModel):
 """Input for weather queries."""
 location: str = Field(description="City name or coordinates")
 units: Literal["celsius", "fahrenheit"] = Field(
 default="celsius",
 description="Temperature unit preference"
 )
 include_forecast: bool = Field(
 default=False,
 description="Include 5-day forecast"
 )

@tool(args_schema=WeatherInput)
def get_weather(location: str, units: str = "celsius", include_forecast: bool = False) -> str:
 """Get current weather and optional forecast."""
 temp = 22 if units == "celsius" else 72
 result = f"Current weather in {location}: {temp} degrees {units[0].upper()}"
 if include_forecast:
 result += "\nNext 5 days: Sunny"
 return result

```

### [​](#reserved-argument-names)Reserved argument names

The following parameter names are reserved and cannot be used as tool arguments. Using these names will cause runtime errors.

| 
 | Parameter name | Purpose
 | `config` | Reserved for passing `RunnableConfig` to tools internally
 | `runtime` | Reserved for `ToolRuntime` parameter (accessing state, context, store)
To access runtime information, use the [`ToolRuntime`](https://reference.langchain.com/python/langchain/tools/#langchain.tools.ToolRuntime) parameter instead of naming your own arguments `config` or `runtime`.

## [​](#access-context)Access context

Tools are most powerful when they can access runtime information like conversation history, user data, and persistent memory. This section covers how to access and update this information from within your tools.
Tools can access runtime information through the [`ToolRuntime`](https://reference.langchain.com/python/langchain/tools/#langchain.tools.ToolRuntime) parameter, which provides:

| 
 | Component | Description | Use case
 | **State** | Short-term memory - mutable data that exists for the current conversation (messages, counters, custom fields) | Access conversation history, track tool call counts
 | **Context** | Immutable configuration passed at invocation time (user IDs, session info) | Personalize responses based on user identity
 | **Store** | Long-term memory - persistent data that survives across conversations | Save user preferences, maintain knowledge base
 | **Stream Writer** | Emit real-time updates during tool execution | Show progress for long-running operations
 | **Execution Info** | Identity and retry information for the current execution (thread ID, run ID, attempt number) | Access thread/run IDs, adjust behavior based on retry state
 | **Server Info** | Server-specific metadata when running on LangGraph Server (assistant ID, graph ID, authenticated user) | Access assistant ID, graph ID, or authenticated user info
 | **Config** | [`RunnableConfig`](https://reference.langchain.com/python/langchain-core/runnables/config/RunnableConfig) for the execution | Access callbacks, tags, and metadata
 | **Tool Call ID** | Unique identifier for the current tool invocation | Correlate tool calls for logs and model invocations

### [​](#short-term-memory-state)Short-term memory (State)

State represents short-term memory that exists for the duration of a conversation. It includes the message history and any custom fields you define in your [graph state](/oss/python/langgraph/graph-api#state).
Add `runtime: ToolRuntime` to your tool signature to access state. This parameter is automatically injected and hidden from the LLM - it won’t appear in the tool’s schema.

#### [​](#access-state)Access state

Tools can access the current conversation state using `runtime.state`:

```
from langchain.tools import tool, ToolRuntime
from langchain.messages import HumanMessage

@tool
def get_last_user_message(runtime: ToolRuntime) -> str:
 """Get the most recent message from the user."""
 messages = runtime.state["messages"]

 # Find the last human message
 for message in reversed(messages):
 if isinstance(message, HumanMessage):
 return message.content

 return "No user messages found"

# Access custom state fields
@tool
def get_user_preference(
 pref_name: str,
 runtime: ToolRuntime
) -> str:
 """Get a user preference value."""
 preferences = runtime.state.get("user_preferences", {})
 return preferences.get(pref_name, "Not set")

```

The `runtime` parameter is hidden from the model. For the example above, the model only sees `pref_name` in the tool schema.

#### [​](#update-state)Update state

Use [`Command`](https://reference.langchain.com/python/langgraph/types/Command) to update the agent’s state. This is useful for tools that need to update custom state fields.
Include a `ToolMessage` in the update so the model can see the result of the tool call:

```
from langchain.agents import AgentState
from langchain.messages import ToolMessage
from langchain.tools import ToolRuntime, tool
from langgraph.types import Command

class CustomState(AgentState):
 user_name: str

@tool
def set_user_name(new_name: str, runtime: ToolRuntime[None, CustomState]) -> Command:
 """Set the user's name in the conversation state."""
 return Command(
 update={
 "user_name": new_name,
 "messages": [
 ToolMessage(
 content=f"User name set to {new_name}.",
 tool_call_id=runtime.tool_call_id,
 )
 ],
 }
 )

```

When tools update state variables, consider defining a [reducer](/oss/python/langgraph/graph-api#reducers) for those fields. Since LLMs can call multiple tools in parallel, a reducer determines how to resolve conflicts when the same state field is updated by concurrent tool calls.

### [​](#context)Context

Context provides immutable configuration data that is passed at invocation time. Use it for user IDs, session details, or application-specific settings that shouldn’t change during a conversation.
While `thread_id` (passed via `config={"configurable": {"thread_id": ...}}`) scopes the *conversation*: message history and checkpoints, `context` carries *per-run* data your tools and middleware read at invocation time. In production you typically pass both together: a stable `thread_id` per conversation, and a `context` object on every invoke.
Access context through `runtime.context`. Pass it alongside a `thread_id` so the conversation is persisted across turns:
GoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from dataclasses import dataclass

from langchain.agents import create_agent
from langchain.tools import tool, ToolRuntime
from langchain_core.utils.uuid import uuid7
from langchain_openai import ChatOpenAI

USER_DATABASE = {
 "user123": {
 "name": "Alice Johnson",
 "account_type": "Premium",
 "balance": 5000,
 "email": "alice@example.com",
 },
 "user456": {
 "name": "Bob Smith",
 "account_type": "Standard",
 "balance": 1200,
 "email": "bob@example.com",
 },
}

@dataclass
class UserContext:
 user_id: str

@tool
def get_account_info(runtime: ToolRuntime[UserContext]) -> str:
 """Get the current user's account information."""
 user_id = runtime.context.user_id

 if user_id in USER_DATABASE:
 user = USER_DATABASE[user_id]
 return (
 f"Account holder: {user['name']}\n"
 f"Type: {user['account_type']}\n"
 f"Balance: ${user['balance']}"
 )
 return "User not found"

model = ChatOpenAI(model="google_genai:gemini-3.5-flash")
agent = create_agent(
 model,
 tools=[get_account_info],
 context_schema=UserContext,
 system_prompt="You are a financial assistant.",
)

result = agent.invoke(
 {"messages": [{"role": "user", "content": "What's my current balance?"}]},
 config={"configurable": {"thread_id": str(uuid7())}},
 context=UserContext(user_id="user123"),
)

```

### [​](#long-term-memory-store)Long-term memory (Store)

The [`BaseStore`](https://reference.langchain.com/python/langchain-core/stores/BaseStore) provides persistent storage that survives across conversations. Unlike state (short-term memory), data saved to the store remains available in future sessions.
Access the store through `runtime.store`. The store uses a namespace/key pattern to organize data:
For production deployments, use a persistent store implementation like [`PostgresStore`](https://reference.langchain.com/python/langgraph/store/#langgraph.store.postgres.PostgresStore) instead of `InMemoryStore`. See the [memory documentation](/oss/python/langgraph/add-memory) for setup details.

```
from typing import Any
from langgraph.store.memory import InMemoryStore
from langchain.agents import create_agent
from langchain.tools import tool, ToolRuntime
from langchain_openai import ChatOpenAI

# Access memory
@tool
def get_user_info(user_id: str, runtime: ToolRuntime) -> str:
 """Look up user info."""
 store = runtime.store
 user_info = store.get(("users",), user_id)
 return str(user_info.value) if user_info else "Unknown user"

# Update memory
@tool
def save_user_info(user_id: str, user_info: dict[str, Any], runtime: ToolRuntime) -> str:
 """Save user info."""
 store = runtime.store
 store.put(("users",), user_id, user_info)
 return "Successfully saved user info."

model = ChatOpenAI(model="gpt-5.4")

store = InMemoryStore()
agent = create_agent(
 model,
 tools=[get_user_info, save_user_info],
 store=store
)

# First session: save user info
agent.invoke({
 "messages": [{"role": "user", "content": "Save the following user: userid: abc123, name: Foo, age: 25, email: foo@langchain.dev"}]
})

# Second session: get user info
agent.invoke({
 "messages": [{"role": "user", "content": "Get user info for user with id 'abc123'"}]
})
# Here is the user info for user with ID "abc123":
# - Name: Foo
# - Age: 25
# - Email: foo@langchain.dev

```
See all 44 lines

### [​](#stream-writer)Stream writer

Stream real-time updates from tools during execution. This is useful for providing progress feedback to users during long-running operations.
Use `runtime.stream_writer` to emit custom updates:

```
from langchain.tools import tool, ToolRuntime

@tool
def get_weather(city: str, runtime: ToolRuntime) -> str:
 """Get weather for a given city."""
 writer = runtime.stream_writer

 # Stream custom updates as the tool executes
 writer(f"Looking up data for city: {city}")
 writer(f"Acquired data for city: {city}")

 return f"It's always sunny in {city}!"

```

If you use `runtime.stream_writer` inside your tool, the tool must be invoked within a LangGraph execution context. See [Streaming](/oss/python/langchain/streaming) for more details.

### [​](#execution-info)Execution info

Access thread ID, run ID, and retry state from within a tool via `runtime.execution_info`:

```
from langchain.tools import tool, ToolRuntime

@tool
def log_execution_context(runtime: ToolRuntime) -> str:
 """Log execution identity information."""
 info = runtime.execution_info
 print(f"Thread: {info.thread_id}, Run: {info.run_id}")
 print(f"Attempt: {info.node_attempt}")
 return "done"

```

Requires `deepagents>=0.5.0` (or `langgraph>=1.1.5`).

### [​](#server-info)Server info

When your tool runs on LangGraph Server, access the assistant ID, graph ID, and authenticated user via `runtime.server_info`:

```
from langchain.tools import tool, ToolRuntime

@tool
def get_assistant_scoped_data(runtime: ToolRuntime) -> str:
 """Fetch data scoped to the current assistant."""
 server = runtime.server_info
 if server is not None:
 print(f"Assistant: {server.assistant_id}, Graph: {server.graph_id}")
 if server.user is not None:
 print(f"User: {server.user.identity}")
 return "done"

```

`server_info` is `None` when the tool is not running on LangGraph Server (e.g., during local development or testing).
Requires `deepagents>=0.5.0` (or `langgraph>=1.1.5`).

## [​](#tool-execution)Tool execution

In LangChain, tools are used by agents (for example via [`create_agent`](https://reference.langchain.com/python/langchain/agents/factory/create_agent)) and tool error handling is configured through [middleware](/oss/python/langchain/middleware).
For LangGraph workflows, tool execution is handled by [`ToolNode`](https://reference.langchain.com/python/langgraph/agents/#langgraph.prebuilt.tool_node.ToolNode). See [ToolNode](/oss/python/langgraph/workflows-agents#toolnode).

### [​](#tool-return-values)Tool return values

You can choose different return values for your tools:

{' ' * (self.list_depth - 1)}- Return a `string` for human-readable results.

{' ' * (self.list_depth - 1)}- Return an `object` for structured results the model should parse.

{' ' * (self.list_depth - 1)}- Return a `Command` with optional message when you need to write to state.

#### [​](#return-a-string)Return a string

Return a string when the tool should provide plain text for the model to read and use in its next response.

```
from langchain.tools import tool

@tool
def get_weather(city: str) -> str:
 """Get weather for a city."""
 return f"It is currently sunny in {city}."

```

Behavior:

{' ' * (self.list_depth - 1)}- The return value is converted to a `ToolMessage`.

{' ' * (self.list_depth - 1)}- The model sees that text and decides what to do next.

{' ' * (self.list_depth - 1)}- No agent state fields are changed unless the model or another tool does so later.

Use this when the result is naturally human-readable text.

#### [​](#return-an-object)Return an object

Return an object (for example, a `dict`) when your tool produces structured data that the model should inspect.

```
from langchain.tools import tool

@tool
def get_weather_data(city: str) -> dict:
 """Get structured weather data for a city."""
 return {
 "city": city,
 "temperature_c": 22,
 "conditions": "sunny",
 }

```

Behavior:

{' ' * (self.list_depth - 1)}- The object is serialized and sent back as tool output.

{' ' * (self.list_depth - 1)}- The model can read specific fields and reason over them.

{' ' * (self.list_depth - 1)}- Like string returns, this does not directly update graph state.

Use this when downstream reasoning benefits from explicit fields instead of free-form text.

#### [​](#return-a-command)Return a Command

Return a [`Command`](https://reference.langchain.com/python/langgraph/types/Command) when the tool needs to update graph state (for example, setting user preferences or app state).
You can return a `Command` with or without including a `ToolMessage`.
If the model needs to see that the tool succeeded (for example, to confirm a preference change), include a `ToolMessage` in the update, using `runtime.tool_call_id` for the `tool_call_id` parameter.

```
from langchain.messages import ToolMessage
from langchain.tools import ToolRuntime, tool
from langgraph.types import Command

@tool
def set_language(language: str, runtime: ToolRuntime) -> Command:
 """Set the preferred response language."""
 return Command(
 update={
 "preferred_language": language,
 "messages": [
 ToolMessage(
 content=f"Language set to {language}.",
 tool_call_id=runtime.tool_call_id,
 )
 ],
 }
 )

```

Behavior:

{' ' * (self.list_depth - 1)}- The command updates state using `update`.

{' ' * (self.list_depth - 1)}- Updated state is available to subsequent steps in the same run.

{' ' * (self.list_depth - 1)}- Use reducers for fields that may be updated by parallel tool calls.

Use this when the tool is not just returning data, but also mutating agent state.

### [​](#error-handling)Error handling

Handle tool errors using LangChain agent [middleware](/oss/python/langchain/middleware) to retry failed tool calls or return custom error messages:
GoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from collections.abc import Callable

from langchain.agents import create_agent
from langchain.agents.middleware import wrap_tool_call
from langchain.messages import ToolMessage
from langchain.tools.tool_node import ToolCallRequest

@wrap_tool_call
def handle_tool_errors(
 request: ToolCallRequest,
 handler: Callable[[ToolCallRequest], ToolMessage],
) -> ToolMessage:
 """Convert tool exceptions into ToolMessages the model can handle."""
 try:
 return handler(request)
 except Exception as e:
 return ToolMessage(
 content=f"Tool error: Please check your input and try again. ({e})",
 tool_call_id=request.tool_call["id"],
 )

agent = create_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[],
 middleware=[handle_tool_errors],
)

```

### [​](#state-injection)State injection

Tools can access the current graph state through [`ToolRuntime`](https://reference.langchain.com/python/langchain/tools/#langchain.tools.ToolRuntime):

```
from langchain.tools import tool, ToolRuntime

@tool
def get_message_count(runtime: ToolRuntime) -> str:
 """Get the number of messages in the conversation."""
 messages = runtime.state["messages"]
 return f"There are {len(messages)} messages."

```

For more details on accessing state, context, and long-term memory from tools, see [Access context](#access-context).

## [​](#dynamic-tool-selection)Dynamic tool selection

With dynamic tools, the set of tools available to the agent is modified at runtime rather than defined all upfront. Not every tool is appropriate for every situation. Too many tools may overwhelm the model (overload context) and increase errors; too few limit capabilities. Dynamic tool selection enables adapting the available toolset based on authentication state, user permissions, feature flags, or conversation stage.
There are two approaches depending on whether tools are known ahead of time:

{' ' * (self.list_depth - 1)}- Filtering pre-registered tools
{' ' * (self.list_depth - 1)}- Runtime tool registrationWhen all possible tools are known at agent creation time, you can pre-register them and dynamically filter which ones are exposed to the model based on state, permissions, or context.
{' ' * (self.list_depth - 1)}- State
{' ' * (self.list_depth - 1)}- Store
{' ' * (self.list_depth - 1)}- Runtime ContextEnable advanced tools only after certain conversation milestones:
```
from langchain.agents import create_agent
from langchain.agents.middleware import wrap_model_call, ModelRequest, ModelResponse
from typing import Callable

@wrap_model_call
def state_based_tools(
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse]
) -> ModelResponse:
 """Filter tools based on conversation State."""
 # Read from State: check if user has authenticated
 state = request.state
 is_authenticated = state.get("authenticated", False)
 message_count = len(state["messages"])

 # Only enable sensitive tools after authentication
 if not is_authenticated:
 tools = [t for t in request.tools if t.name.startswith("public_")]
 request = request.override(tools=tools)
 elif message_count < 5:
 # Limit tools early in conversation
 tools = [t for t in request.tools if t.name != "advanced_search"]
 request = request.override(tools=tools)

 return handler(request)

agent = create_agent(
 model="gpt-5.4",
 tools=[public_search, private_search, advanced_search],
 middleware=[state_based_tools]
)

```
Filter tools based on user preferences or feature flags in Store:
```
from dataclasses import dataclass
from langchain.agents import create_agent
from langchain.agents.middleware import wrap_model_call, ModelRequest, ModelResponse
from typing import Callable
from langgraph.store.memory import InMemoryStore

@dataclass
class Context:
 user_id: str

@wrap_model_call
def store_based_tools(
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse]
) -> ModelResponse:
 """Filter tools based on Store preferences."""
 user_id = request.runtime.context.user_id

 # Read from Store: get user's enabled features
 store = request.runtime.store
 feature_flags = store.get(("features",), user_id)

 if feature_flags:
 enabled_features = feature_flags.value.get("enabled_tools", [])
 # Only include tools that are enabled for this user
 tools = [t for t in request.tools if t.name in enabled_features]
 request = request.override(tools=tools)

 return handler(request)

agent = create_agent(
 model="gpt-5.4",
 tools=[search_tool, analysis_tool, export_tool],
 middleware=[store_based_tools],
 context_schema=Context,
 store=InMemoryStore()
)

```
Filter tools based on user permissions from Runtime Context:
```
from dataclasses import dataclass
from langchain.agents import create_agent
from langchain.agents.middleware import wrap_model_call, ModelRequest, ModelResponse
from typing import Callable

@dataclass
class Context:
 user_role: str

@wrap_model_call
def context_based_tools(
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse]
) -> ModelResponse:
 """Filter tools based on Runtime Context permissions."""
 # Read from Runtime Context: get user role
 if request.runtime is None or request.runtime.context is None:
 # If no context provided, default to viewer (most restrictive)
 user_role = "viewer"
 else:
 user_role = request.runtime.context.user_role

 if user_role == "admin":
 # Admins get all tools
 pass
 elif user_role == "editor":
 # Editors can't delete
 tools = [t for t in request.tools if t.name != "delete_data"]
 request = request.override(tools=tools)
 else:
 # Viewers get read-only tools
 tools = [t for t in request.tools if t.name.startswith("read_")]
 request = request.override(tools=tools)

 return handler(request)

agent = create_agent(
 model="gpt-5.4",
 tools=[read_data, write_data, delete_data],
 middleware=[context_based_tools],
 context_schema=Context
)

```
This approach is best when:

{' ' * (self.list_depth - 1)}- All possible tools are known at compile/startup time

{' ' * (self.list_depth - 1)}- You want to filter based on permissions, feature flags, or conversation state

{' ' * (self.list_depth - 1)}- Tools are static but their availability is dynamic
See [Dynamically selecting tools](/oss/python/langchain/middleware/custom#dynamically-selecting-tools) for more examples.When tools are discovered or created at runtime (e.g., loaded from an MCP server, generated based on user data, or fetched from a remote registry), you need to both register the tools and handle their execution dynamically.This requires two middleware hooks:

{' ' * (self.list_depth - 1)}- `wrap_model_call` - Add the dynamic tools to the request

{' ' * (self.list_depth - 1)}- `wrap_tool_call` - Handle execution of the dynamically added tools

```
from langchain.tools import tool
from langchain.agents import create_agent
from langchain.agents.middleware import AgentMiddleware, ModelRequest, ToolCallRequest

# A tool that will be added dynamically at runtime
@tool
def calculate_tip(bill_amount: float, tip_percentage: float = 20.0) -> str:
 """Calculate the tip amount for a bill."""
 tip = bill_amount * (tip_percentage / 100)
 return f"Tip: ${tip:.2f}, Total: ${bill_amount + tip:.2f}"

class DynamicToolMiddleware(AgentMiddleware):
 """Middleware that registers and handles dynamic tools."""

 def wrap_model_call(self, request: ModelRequest, handler):
 # Add dynamic tool to the request
 # This could be loaded from an MCP server, database, etc.
 updated = request.override(tools=[*request.tools, calculate_tip])
 return handler(updated)

 def wrap_tool_call(self, request: ToolCallRequest, handler):
 # Handle execution of the dynamic tool
 if request.tool_call["name"] == "calculate_tip":
 return handler(request.override(tool=calculate_tip))
 return handler(request)

agent = create_agent(
 model="gpt-4o",
 tools=[get_weather], # Only static tools registered here
 middleware=[DynamicToolMiddleware()],
)

# The agent can now use both get_weather AND calculate_tip
result = agent.invoke({
 "messages": [{"role": "user", "content": "Calculate a 20% tip on $85"}]
})

```
This approach is best when:

{' ' * (self.list_depth - 1)}- Tools are discovered at runtime (e.g., from an MCP server)

{' ' * (self.list_depth - 1)}- Tools are generated dynamically based on user data or configuration

{' ' * (self.list_depth - 1)}- You’re integrating with external tool registries
The `wrap_tool_call` hook is required for runtime-registered tools because the agent needs to know how to execute tools that weren’t in the original tool list. Without it, the agent won’t know how to invoke the dynamically added tool.

## [​](#prebuilt-tools)Prebuilt tools

LangChain provides a large collection of prebuilt tools and toolkits for common tasks like web search, code interpretation, database access, and more. These ready-to-use tools can be directly integrated into your agents without writing custom code.
See the [tools and toolkits](/oss/python/integrations/tools) integration page for a complete list of available tools organized by category.

## [​](#server-side-tool-use)Server-side tool use

Some chat models feature built-in tools that are executed server-side by the model provider. These include capabilities like web search and code interpreters that don’t require you to define or host the tool logic.
Refer to the individual [chat model integration pages](/oss/python/integrations/providers) and the [tool calling documentation](/oss/python/langchain/models#server-side-tool-use) for details on enabling and using these built-in tools.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/tools.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[MessagesPrevious](/oss/python/langchain/messages)[Short-term memoryNext](/oss/python/langchain/short-term-memory)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
