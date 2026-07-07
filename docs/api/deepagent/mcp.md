# Model Context Protocol (MCP)

[Advanced usage](/oss/python/langchain/guardrails)

# Model Context Protocol (MCP)
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.[Model Context Protocol (MCP)](https://modelcontextprotocol.io/introduction) is an open protocol that standardizes how applications provide tools and context to LLMs. LangChain agents can use tools defined on MCP servers using the [`langchain-mcp-adapters`](https://github.com/langchain-ai/langchain-mcp-adapters) library.

## [​](#quickstart)Quickstart

Install the `langchain-mcp-adapters` library:
pipuv
```
pip install langchain-mcp-adapters

```

`langchain-mcp-adapters` enables agents to use tools defined across one or more MCP servers.
`MultiServerMCPClient` is **stateless by default**. Each tool invocation creates a fresh MCP `ClientSession`, executes the tool, and then cleans up. See the [stateful sessions](#stateful-sessions) section for more details.
Accessing multiple MCP servers
```
import asyncio
from langchain_mcp_adapters.client import MultiServerMCPClient 
from langchain.agents import create_agent

async def main():
 client = MultiServerMCPClient(
 {
 "math": {
 "transport": "stdio", # Local subprocess communication
 "command": "python",
 # Absolute path to your math_server.py file
 "args": ["/path/to/math_server.py"],
 },
 "weather": {
 "transport": "http", # HTTP-based remote server
 # Ensure you start your weather server on port 8000
 "url": "http://localhost:8000/mcp",
 }
 }
 )

 tools = await client.get_tools()
 agent = create_agent(
 "claude-sonnet-4-6",
 tools 
 )
 math_response = await agent.ainvoke(
 {"messages": [{"role": "user", "content": "what's (3 + 5) x 12?"}]}
 )
 weather_response = await agent.ainvoke(
 {"messages": [{"role": "user", "content": "what is the weather in nyc?"}]}
 )
 print(math_response)
 print(weather_response)

if __name__ == "__main__":
 asyncio.run(main())

```

Trace MCP tool calls alongside your agent’s reasoning steps with [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-mcp). Follow the [tracing quickstart](/langsmith/trace-with-langchain) to get set up.

## [​](#custom-servers)Custom servers

To create a custom MCP server, use the [FastMCP](https://gofastmcp.com/getting-started/welcome) library:
pipuv
```
pip install fastmcp

```

To test your agent with MCP tool servers, use the following examples:
Math server (stdio transport)Weather server (streamable HTTP transport)
```
from fastmcp import FastMCP

mcp = FastMCP("Math")

@mcp.tool()
def add(a: int, b: int) -> int:
 """Add two numbers"""
 return a + b

@mcp.tool()
def multiply(a: int, b: int) -> int:
 """Multiply two numbers"""
 return a * b

if __name__ == "__main__":
 mcp.run(transport="stdio")

```

## [​](#transports)Transports

MCP supports different transport mechanisms for client-server communication.

### [​](#http)HTTP

The `http` transport (also referred to as `streamable-http`) uses HTTP requests for client-server communication. See the [MCP HTTP transport specification](https://modelcontextprotocol.io/specification/2025-03-26/basic/transports#streamable-http) for more details.

```
client = MultiServerMCPClient(
 {
 "weather": {
 "transport": "http",
 "url": "http://localhost:8000/mcp",
 }
 }
)

```

#### [​](#passing-headers)Passing headers

When connecting to MCP servers over HTTP, you can include custom headers (e.g., for authentication or tracing) using the `headers` field in the connection configuration. This is supported for `sse` (deprecated by MCP spec) and `streamable_http` transports.
Passing headers with MultiServerMCPClient
```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain.agents import create_agent

client = MultiServerMCPClient(
 {
 "weather": {
 "transport": "http",
 "url": "http://localhost:8000/mcp",
 "headers": {
 "Authorization": "Bearer YOUR_TOKEN",
 "X-Custom-Header": "custom-value"
 },
 }
 }
)
tools = await client.get_tools()
agent = create_agent("openai:gpt-5.4", tools)
response = await agent.ainvoke({"messages": "what is the weather in nyc?"})

```

#### [​](#authentication)Authentication

The `langchain-mcp-adapters` library uses the official [MCP SDK](https://github.com/modelcontextprotocol/python-sdk) under the hood, which allows you to provide a custom authentication mechanism by implementing the `httpx.Auth` interface.

```
from langchain_mcp_adapters.client import MultiServerMCPClient

client = MultiServerMCPClient(
 {
 "weather": {
 "transport": "http",
 "url": "http://localhost:8000/mcp",
 "auth": auth,
 }
 }
)

```

{' ' * (self.list_depth - 1)}- [Example custom auth implementation](https://github.com/modelcontextprotocol/python-sdk/blob/main/examples/clients/simple-auth-client/mcp_simple_auth_client/main.py)

{' ' * (self.list_depth - 1)}- [Built-in OAuth flow](https://github.com/modelcontextprotocol/python-sdk/blob/main/src/mcp/client/auth/oauth2.py#L216)

### [​](#stdio)stdio

Client launches server as a subprocess and communicates via standard input/output. Best for local tools and simple setups.
Unlike HTTP transports, `stdio` connections are inherently **stateful**: the subprocess persists for the lifetime of the client connection. However, when using `MultiServerMCPClient` without explicit session management, each tool call still creates a new session. See [stateful sessions](#stateful-sessions) for managing persistent connections.

```
client = MultiServerMCPClient(
 {
 "math": {
 "transport": "stdio",
 "command": "python",
 "args": ["/path/to/math_server.py"],
 }
 }
)

```

## [​](#stateful-sessions)Stateful sessions

By default, `MultiServerMCPClient` is **stateless**: each tool invocation creates a fresh MCP session, executes the tool, and then cleans up.
If you need to control the [lifecycle](https://modelcontextprotocol.io/specification/2025-03-26/basic/lifecycle) of an MCP session (for example, when working with a stateful server that maintains context across tool calls), you can create a persistent `ClientSession` using `client.session()`.
Using MCP ClientSession for stateful tool usage
```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.tools import load_mcp_tools
from langchain.agents import create_agent

client = MultiServerMCPClient({...})

# Create a session explicitly
async with client.session("server_name") as session:
 # Pass the session to load tools, resources, or prompts
 tools = await load_mcp_tools(session)
 agent = create_agent(
 "google_genai:gemini-3.5-flash",
 tools
 )

```

## [​](#core-features)Core features

### [​](#tools)Tools

[Tools](https://modelcontextprotocol.io/docs/concepts/tools) allow MCP servers to expose executable functions that LLMs can invoke to perform actions—such as querying databases, calling APIs, or interacting with external systems. LangChain converts MCP tools into LangChain [tools](/oss/python/langchain/tools), making them directly usable in any LangChain agent or workflow.

#### [​](#loading-tools)Loading tools

Use `client.get_tools()` to retrieve tools from MCP servers and pass them to your agent:

```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain.agents import create_agent

client = MultiServerMCPClient({...})
tools = await client.get_tools()
agent = create_agent("claude-sonnet-4-6", tools)

```

#### [​](#structured-content)Structured content

MCP tools can return [structured content](https://modelcontextprotocol.io/specification/2025-03-26/server/tools#structured-content) alongside the human-readable text response. This is useful when a tool needs to return machine-parseable data (like JSON) in addition to text that gets shown to the model.
When an MCP tool returns `structuredContent`, the adapter wraps it in an [`MCPToolArtifact`](https://reference.langchain.com/python/langchain_mcp_adapters/#langchain_mcp_adapters.tools.MCPToolArtifact) and returns it as the tool’s artifact. You can access this using the `artifact` field on the `ToolMessage`. You can also use [interceptors](#tool-interceptors) to process or transform structured content automatically.
**Extracting structured content from artifact**
After invoking your agent, you can access the structured content from tool messages in the response:

```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain.agents import create_agent
from langchain.messages import ToolMessage

client = MultiServerMCPClient({...})
tools = await client.get_tools()
agent = create_agent("claude-sonnet-4-6", tools)

result = await agent.ainvoke(
 {"messages": [{"role": "user", "content": "Get data from the server"}]}
)

# Extract structured content from tool messages
for message in result["messages"]:
 if isinstance(message, ToolMessage) and message.artifact:
 structured_content = message.artifact["structured_content"]

```

**Appending structured content via interceptor**
If you want structured content to be visible in the conversation history (visible to the model), you can use an [interceptor](#tool-interceptors) to automatically append structured content to the tool result:

```
import json

from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.interceptors import MCPToolCallRequest
from mcp.types import TextContent

async def append_structured_content(request: MCPToolCallRequest, handler):
 """Append structured content from artifact to tool message."""
 result = await handler(request)
 if result.structuredContent:
 result.content += [
 TextContent(type="text", text=json.dumps(result.structuredContent)),
 ]
 return result

client = MultiServerMCPClient({...}, tool_interceptors=[append_structured_content])

```

#### [​](#multimodal-tool-content)Multimodal tool content

MCP tools can return [multimodal content](https://modelcontextprotocol.io/specification/2025-03-26/server/tools#tool-result) (images, text, etc.) in their responses. When an MCP server returns content with multiple parts (e.g., text and images), the adapter converts them to LangChain’s [standard content blocks](/oss/python/langchain/messages#standard-content-blocks). You can access the standardized representation via the `content_blocks` property on the `ToolMessage`:

```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain.agents import create_agent

client = MultiServerMCPClient({...})
tools = await client.get_tools()
agent = create_agent("claude-sonnet-4-6", tools)

result = await agent.ainvoke(
 {"messages": [{"role": "user", "content": "Take a screenshot of the current page"}]}
)

# Access multimodal content from tool messages
for message in result["messages"]:
 if message.type == "tool":
 # Raw content in provider-native format
 print(f"Raw content: {message.content}")

 # Standardized content blocks #
 for block in message.content_blocks:
 if block["type"] == "text":
 print(f"Text: {block['text']}")
 elif block["type"] == "image":
 print(f"Image URL: {block.get('url')}")
 print(f"Image base64: {block.get('base64', '')[:50]}...")

```

This allows you to handle multimodal tool responses in a provider-agnostic way, regardless of how the underlying MCP server formats its content.

### [​](#resources)Resources

[Resources](https://modelcontextprotocol.io/docs/concepts/resources) allow MCP servers to expose data—such as files, database records, or API responses—that can be read by clients. LangChain converts MCP resources into [Blob](https://reference.langchain.com/python/langchain_core/documents/#langchain_core.documents.base.Blob) objects, which provide a unified interface for handling both text and binary content.

#### [​](#loading-resources)Loading resources

Use `client.get_resources()` to load resources from an MCP server:

```
from langchain_mcp_adapters.client import MultiServerMCPClient

client = MultiServerMCPClient({...})

# Load all resources from a server
blobs = await client.get_resources("server_name")

# Or load specific resources by URI
blobs = await client.get_resources("server_name", uris=["file:///path/to/file.txt"])

for blob in blobs:
 print(f"URI: {blob.metadata['uri']}, MIME type: {blob.mimetype}")
 print(blob.as_string()) # For text content

```

You can also use [`load_mcp_resources`](https://reference.langchain.com/python/langchain_mcp_adapters/#langchain_mcp_adapters.resources.load_mcp_resources) directly with a session for more control:

```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.resources import load_mcp_resources

client = MultiServerMCPClient({...})

async with client.session("server_name") as session:
 # Load all resources
 blobs = await load_mcp_resources(session)

 # Or load specific resources by URI
 blobs = await load_mcp_resources(session, uris=["file:///path/to/file.txt"])

```

### [​](#prompts)Prompts

[Prompts](https://modelcontextprotocol.io/docs/concepts/prompts) allow MCP servers to expose reusable prompt templates that can be retrieved and used by clients. LangChain converts MCP prompts into [messages](/oss/python/langchain/messages), making them easy to integrate into chat-based workflows.

#### [​](#loading-prompts)Loading prompts

Use `client.get_prompt()` to load a prompt from an MCP server:

```
from langchain_mcp_adapters.client import MultiServerMCPClient

client = MultiServerMCPClient({...})

# Load a prompt by name
messages = await client.get_prompt("server_name", "summarize")

# Load a prompt with arguments
messages = await client.get_prompt(
 "server_name",
 "code_review",
 arguments={"language": "python", "focus": "security"}
)

# Use the messages in your workflow
for message in messages:
 print(f"{message.type}: {message.content}")

```

You can also use [`load_mcp_prompt`](https://reference.langchain.com/python/langchain_mcp_adapters/#langchain_mcp_adapters.prompts.load_mcp_prompt) directly with a session for more control:

```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.prompts import load_mcp_prompt

client = MultiServerMCPClient({...})

async with client.session("server_name") as session:
 # Load a prompt by name
 messages = await load_mcp_prompt(session, "summarize")

 # Load a prompt with arguments
 messages = await load_mcp_prompt(
 session,
 "code_review",
 arguments={"language": "python", "focus": "security"}
 )

```

## [​](#advanced-features)Advanced features

### [​](#tool-interceptors)Tool interceptors

MCP servers run as separate processes—they can’t access LangGraph runtime information like the [store](/oss/python/langgraph/persistence#memory-store), [context](/oss/python/langchain/context-engineering), or agent state. **Interceptors bridge this gap** by giving you access to this runtime context during MCP tool execution.
Interceptors also provide middleware-like control over tool calls: you can modify requests, implement retries, add headers dynamically, or short-circuit execution entirely.

| 
 | Section | Description
 | [Accessing runtime context](#accessing-runtime-context) | Read user IDs, API keys, store data, and agent state
 | [State updates and commands](#state-updates-and-commands) | Update agent state or control graph flow with `Command`
 | [Writing interceptors](#custom-interceptors) | Patterns for modifying requests, composing interceptors, and error handling

#### [​](#accessing-runtime-context)Accessing runtime context

When MCP tools are used within a LangChain agent (via `create_agent`), interceptors receive access to the `ToolRuntime` context. This provides access to the tool call ID, state, config, and store—enabling powerful patterns for accessing user data, persisting information, and controlling agent behavior.

{' ' * (self.list_depth - 1)}- Runtime context
{' ' * (self.list_depth - 1)}- Store
{' ' * (self.list_depth - 1)}- State
{' ' * (self.list_depth - 1)}- Tool call IDAccess user-specific configuration like user IDs, API keys, or permissions that are passed at invocation time:Inject user context into MCP tool calls
```
from dataclasses import dataclass
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.interceptors import MCPToolCallRequest
from langchain.agents import create_agent

@dataclass
class Context:
 user_id: str
 api_key: str

async def inject_user_context(
 request: MCPToolCallRequest,
 handler,
):
 """Inject user credentials into MCP tool calls."""
 runtime = request.runtime
 user_id = runtime.context.user_id 
 api_key = runtime.context.api_key 

 # Add user context to tool arguments
 modified_request = request.override(
 args={**request.args, "user_id": user_id}
 )
 return await handler(modified_request)

client = MultiServerMCPClient(
 {...},
 tool_interceptors=[inject_user_context],
)
tools = await client.get_tools()
agent = create_agent("gpt-5.4", tools, context_schema=Context)

# Invoke with user context
result = await agent.ainvoke(
 {"messages": [{"role": "user", "content": "Search my orders"}]},
 context={"user_id": "user_123", "api_key": "sk-..."}
)

```
Access long-term memory to retrieve user preferences or persist data across conversations:Access user preferences from store
```
from dataclasses import dataclass
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.interceptors import MCPToolCallRequest
from langchain.agents import create_agent
from langgraph.store.memory import InMemoryStore

@dataclass
class Context:
 user_id: str

async def personalize_search(
 request: MCPToolCallRequest,
 handler,
):
 """Personalize MCP tool calls using stored preferences."""
 runtime = request.runtime
 user_id = runtime.context.user_id
 store = runtime.store 

 # Read user preferences from store
 prefs = store.get(("preferences",), user_id)

 if prefs and request.name == "search":
 # Apply user's preferred language and result limit
 modified_args = {
 **request.args,
 "language": prefs.value.get("language", "en"),
 "limit": prefs.value.get("result_limit", 10),
 }
 request = request.override(args=modified_args)

 return await handler(request)

client = MultiServerMCPClient(
 {...},
 tool_interceptors=[personalize_search],
)
tools = await client.get_tools()
agent = create_agent(
 "gpt-5.4",
 tools,
 context_schema=Context,
 store=InMemoryStore()
)

```
Access conversation state to make decisions based on the current session:Filter tools based on authentication state
```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.interceptors import MCPToolCallRequest
from langchain.messages import ToolMessage

async def require_authentication(
 request: MCPToolCallRequest,
 handler,
):
 """Block sensitive MCP tools if user is not authenticated."""
 runtime = request.runtime
 state = runtime.state 
 is_authenticated = state.get("authenticated", False)

 sensitive_tools = ["delete_file", "update_settings", "export_data"]

 if request.name in sensitive_tools and not is_authenticated:
 # Return error instead of calling tool
 return ToolMessage(
 content="Authentication required. Please log in first.",
 tool_call_id=runtime.tool_call_id,
 )

 return await handler(request)

client = MultiServerMCPClient(
 {...},
 tool_interceptors=[require_authentication],
)

```
Access the tool call ID to return properly formatted responses or track tool executions:Return custom responses with tool call ID
```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.interceptors import MCPToolCallRequest
from langchain.messages import ToolMessage

async def rate_limit_interceptor(
 request: MCPToolCallRequest,
 handler,
):
 """Rate limit expensive MCP tool calls."""
 runtime = request.runtime
 tool_call_id = runtime.tool_call_id 

 # Check rate limit (simplified example)
 if is_rate_limited(request.name):
 return ToolMessage(
 content="Rate limit exceeded. Please try again later.",
 tool_call_id=tool_call_id,
 )

 result = await handler(request)

 # Log successful tool call
 log_tool_execution(tool_call_id, request.name, success=True)

 return result

client = MultiServerMCPClient(
 {...},
 tool_interceptors=[rate_limit_interceptor],
)

```

For more context engineering patterns, see [Context engineering](/oss/python/langchain/context-engineering) and [Tools](/oss/python/langchain/tools).

#### [​](#state-updates-and-commands)State updates and commands

Interceptors can return `Command` objects to update agent state or control graph execution flow. This is useful for tracking task progress, switching between agents, or ending execution early.
Mark task complete and switch agents
```
from langchain.agents import AgentState, create_agent
from langchain_mcp_adapters.interceptors import MCPToolCallRequest
from langchain.messages import ToolMessage
from langgraph.types import Command

async def handle_task_completion(
 request: MCPToolCallRequest,
 handler,
):
 """Mark task complete and hand off to summary agent."""
 result = await handler(request)

 if request.name == "submit_order":
 return Command(
 update={
 "messages": [result] if isinstance(result, ToolMessage) else [],
 "task_status": "completed",
 },
 goto="summary_agent",
 )

 return result

```

Use `Command` with `goto="__end__"` to end execution early:
End agent run on completion
```
async def end_on_success(
 request: MCPToolCallRequest,
 handler,
):
 """End agent run when task is marked complete."""
 result = await handler(request)

 if request.name == "mark_complete":
 return Command(
 update={"messages": [result], "status": "done"},
 goto="__end__",
 )

 return result

```

#### [​](#custom-interceptors)Custom interceptors

Interceptors are async functions that wrap tool execution, enabling request/response modification, retry logic, and other cross-cutting concerns. They follow an “onion” pattern where the first interceptor in the list is the outermost layer.
**Basic pattern**
An interceptor is an async function that receives a request and a handler. You can modify the request before calling the handler, modify the response after, or skip the handler entirely.
Basic interceptor pattern
```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.interceptors import MCPToolCallRequest

async def logging_interceptor(
 request: MCPToolCallRequest,
 handler,
):
 """Log tool calls before and after execution."""
 print(f"Calling tool: {request.name} with args: {request.args}")
 result = await handler(request)
 print(f"Tool {request.name} returned: {result}")
 return result

client = MultiServerMCPClient(
 {"math": {"transport": "stdio", "command": "python", "args": ["/path/to/server.py"]}},
 tool_interceptors=[logging_interceptor],
)

```

**Modifying requests**
Use `request.override()` to create a modified request. This follows an immutable pattern, leaving the original request unchanged.
Modifying tool arguments
```
async def double_args_interceptor(
 request: MCPToolCallRequest,
 handler,
):
 """Double all numeric arguments before execution."""
 modified_args = {k: v * 2 for k, v in request.args.items()}
 modified_request = request.override(args=modified_args)
 return await handler(modified_request)

# Original call: add(a=2, b=3) becomes add(a=4, b=6)

```

**Modifying headers at runtime**
Interceptors can modify HTTP headers dynamically based on the request context:
Dynamic header modification
```
async def auth_header_interceptor(
 request: MCPToolCallRequest,
 handler,
):
 """Add authentication headers based on the tool being called."""
 token = get_token_for_tool(request.name)
 modified_request = request.override(
 headers={"Authorization": f"Bearer {token}"}
 )
 return await handler(modified_request)

```

**Composing interceptors**
Multiple interceptors compose in “onion” order—the first interceptor in the list is the outermost layer:
Composing multiple interceptors
```
async def outer_interceptor(request, handler):
 print("outer: before")
 result = await handler(request)
 print("outer: after")
 return result

async def inner_interceptor(request, handler):
 print("inner: before")
 result = await handler(request)
 print("inner: after")
 return result

client = MultiServerMCPClient(
 {...},
 tool_interceptors=[outer_interceptor, inner_interceptor],
)

# Execution order:
# outer: before -> inner: before -> tool execution -> inner: after -> outer: after

```

**Error handling**
Use interceptors to catch tool execution errors and implement retry logic:
Retry on error
```
import asyncio

async def retry_interceptor(
 request: MCPToolCallRequest,
 handler,
 max_retries: int = 3,
 delay: float = 1.0,
):
 """Retry failed tool calls with exponential backoff."""
 last_error = None
 for attempt in range(max_retries):
 try:
 return await handler(request)
 except Exception as e:
 last_error = e
 if attempt < max_retries - 1:
 wait_time = delay * (2 ** attempt) # Exponential backoff
 print(f"Tool {request.name} failed (attempt {attempt + 1}), retrying in {wait_time}s...")
 await asyncio.sleep(wait_time)
 raise last_error

client = MultiServerMCPClient(
 {...},
 tool_interceptors=[retry_interceptor],
)

```

You can also catch specific error types and return fallback values:
Error handling with fallback
```
async def fallback_interceptor(
 request: MCPToolCallRequest,
 handler,
):
 """Return a fallback value if tool execution fails."""
 try:
 return await handler(request)
 except TimeoutError:
 return f"Tool {request.name} timed out. Please try again later."
 except ConnectionError:
 return f"Could not connect to {request.name} service. Using cached data."

```

### [​](#progress-notifications)Progress notifications

Subscribe to progress updates for long-running tool executions:
Progress callback
```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.callbacks import Callbacks, CallbackContext

async def on_progress(
 progress: float,
 total: float | None,
 message: str | None,
 context: CallbackContext,
):
 """Handle progress updates from MCP servers."""
 percent = (progress / total * 100) if total else progress
 tool_info = f" ({context.tool_name})" if context.tool_name else ""
 print(f"[{context.server_name}{tool_info}] Progress: {percent:.1f}% - {message}")

client = MultiServerMCPClient(
 {...},
 callbacks=Callbacks(on_progress=on_progress),
)

```

The `CallbackContext` provides:

{' ' * (self.list_depth - 1)}- `server_name`: Name of the MCP server

{' ' * (self.list_depth - 1)}- `tool_name`: Name of the tool being executed (available during tool calls)

### [​](#logging)Logging

The MCP protocol supports [logging](https://modelcontextprotocol.io/specification/2025-03-26/server/utilities/logging#log-levels) notifications from servers. Use the `Callbacks` class to subscribe to these events.
Logging callback
```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.callbacks import Callbacks, CallbackContext
from mcp.types import LoggingMessageNotificationParams

async def on_logging_message(
 params: LoggingMessageNotificationParams,
 context: CallbackContext,
):
 """Handle log messages from MCP servers."""
 print(f"[{context.server_name}] {params.level}: {params.data}")

client = MultiServerMCPClient(
 {...},
 callbacks=Callbacks(on_logging_message=on_logging_message),
)

```

### [​](#elicitation)Elicitation

[Elicitation](https://modelcontextprotocol.io/specification/2025-11-25/client/elicitation#elicitation) allows MCP servers to request additional input from users during tool execution. Instead of requiring all inputs upfront, servers can interactively ask for information as needed.

#### [​](#server-setup)Server setup

Define a tool that uses `ctx.elicit()` to request user input with a schema:
MCP server with elicitation
```
from pydantic import BaseModel
from mcp.server.fastmcp import Context, FastMCP

server = FastMCP("Profile")

class UserDetails(BaseModel):
 email: str
 age: int

@server.tool()
async def create_profile(name: str, ctx: Context) -> str:
 """Create a user profile, requesting details via elicitation."""
 result = await ctx.elicit(
 message=f"Please provide details for {name}'s profile:",
 schema=UserDetails,
 )
 if result.action == "accept" and result.data:
 return f"Created profile for {name}: email={result.data.email}, age={result.data.age}"
 if result.action == "decline":
 return f"User declined. Created minimal profile for {name}."
 return "Profile creation cancelled."

if __name__ == "__main__":
 server.run(transport="http")

```

#### [​](#client-setup)Client setup

Handle elicitation requests by providing a callback to `MultiServerMCPClient`:
Handling elicitation requests
```
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_mcp_adapters.callbacks import Callbacks, CallbackContext
from mcp.shared.context import RequestContext
from mcp.types import ElicitRequestParams, ElicitResult

async def on_elicitation(
 mcp_context: RequestContext,
 params: ElicitRequestParams,
 context: CallbackContext,
) -> ElicitResult:
 """Handle elicitation requests from MCP servers."""
 # In a real application, you would prompt the user for input
 # based on params.message and params.requestedSchema
 return ElicitResult(
 action="accept",
 content={"email": "user@example.com", "age": 25},
 )

client = MultiServerMCPClient(
 {
 "profile": {
 "url": "http://localhost:8000/mcp",
 "transport": "http",
 }
 },
 callbacks=Callbacks(on_elicitation=on_elicitation),
)

```

#### [​](#response-actions)Response actions

The elicitation callback can return one of three actions:

| 
 | Action | Description
 | `accept` | User provided valid input. Include the data in the `content` field.
 | `decline` | User chose not to provide the requested information.
 | `cancel` | User cancelled the operation entirely.
Response action examples
```
# Accept with data
ElicitResult(action="accept", content={"email": "user@example.com", "age": 25})

# Decline (user doesn't want to provide info)
ElicitResult(action="decline")

# Cancel (abort the operation)
ElicitResult(action="cancel")

```

## [​](#additional-resources)Additional resources

{' ' * (self.list_depth - 1)}- [MCP documentation](https://modelcontextprotocol.io/introduction)

{' ' * (self.list_depth - 1)}- [MCP Transport documentation](https://modelcontextprotocol.io/docs/concepts/transports)

{' ' * (self.list_depth - 1)}- [`langchain-mcp-adapters`](https://github.com/langchain-ai/langchain-mcp-adapters)

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/mcp.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Context engineering in agentsPrevious](/oss/python/langchain/context-engineering)[Human-in-the-loopNext](/oss/python/langchain/human-in-the-loop)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
