# Custom middleware

[Middleware](/oss/python/langchain/middleware/overview)

# Custom middleware
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Build custom middleware by implementing hooks that run at specific points in the agent execution flow.

## [​](#hooks)Hooks

Middleware provides two styles of hooks to intercept agent execution:

## Node-style hooks
Run sequentially at specific execution points.

## Wrap-style hooks
Run around each model or tool call.

### [​](#node-style-hooks)Node-style hooks

Run sequentially at specific execution points. Use for logging, validation, and state updates.
Choose the hooks your middleware needs. You can choose between node-style hooks and wrap-style hooks.
**Node-style hooks** run at specific execution points:

| 
 | Hook | When it runs
 | `before_agent` | Before agent starts (once per invocation)
 | `before_model` | Before each model call
 | `after_model` | After each model response
 | `after_agent` | After agent completes (once per invocation)
**Wrap-style hooks** run around each call, giving you control over execution:

| 
 | Hook | When it runs
 | `wrap_model_call` | Around each model call
 | `wrap_tool_call` | Around each tool call
**Example:**

{' ' * (self.list_depth - 1)}- Decorator
{' ' * (self.list_depth - 1)}- Class
```
from langchain.agents.middleware import before_model, after_model, AgentState
from langchain.messages import AIMessage
from langgraph.runtime import Runtime
from typing import Any

@before_model(can_jump_to=["end"])
def check_message_limit(state: AgentState, runtime: Runtime) -> dict[str, Any] | None:
 if len(state["messages"]) >= 50:
 return {
 "messages": [AIMessage("Conversation limit reached.")],
 "jump_to": "end"
 }
 return None

@after_model
def log_response(state: AgentState, runtime: Runtime) -> dict[str, Any] | None:
 print(f"Model returned: {state['messages'][-1].content}")
 return None

```

```
from langchain.agents.middleware import AgentMiddleware, AgentState, hook_config
from langchain.messages import AIMessage
from langgraph.runtime import Runtime
from typing import Any

class MessageLimitMiddleware(AgentMiddleware):
 def __init__(self, max_messages: int = 50):
 super().__init__()
 self.max_messages = max_messages

 @hook_config(can_jump_to=["end"])
 def before_model(self, state: AgentState, runtime: Runtime) -> dict[str, Any] | None:
 if len(state["messages"]) >= self.max_messages:
 return {
 "messages": [AIMessage("Conversation limit reached.")],
 "jump_to": "end"
 }
 return None

 def after_model(self, state: AgentState, runtime: Runtime) -> dict[str, Any] | None:
 print(f"Model returned: {state['messages'][-1].content}")
 return None

```

### [​](#wrap-style-hooks)Wrap-style hooks

Intercept execution and control when the handler is called. Use for retries, caching, and transformation.
You decide if the handler is called zero times (short-circuit), once (normal flow), or multiple times (retry logic).
**Available hooks:**

{' ' * (self.list_depth - 1)}- `wrap_model_call` - Around each model call

{' ' * (self.list_depth - 1)}- `wrap_tool_call` - Around each tool call

**Example:**

{' ' * (self.list_depth - 1)}- Decorator
{' ' * (self.list_depth - 1)}- Class
```
from langchain.agents.middleware import wrap_model_call, ModelRequest, ModelResponse
from typing import Callable

@wrap_model_call
def retry_model(
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
) -> ModelResponse:
 for attempt in range(3):
 try:
 return handler(request)
 except Exception as e:
 if attempt == 2:
 raise
 print(f"Retry {attempt + 1}/3 after error: {e}")

```

```
from langchain.agents.middleware import AgentMiddleware, ModelRequest, ModelResponse
from typing import Callable

class RetryMiddleware(AgentMiddleware):
 def __init__(self, max_retries: int = 3):
 super().__init__()
 self.max_retries = max_retries

 def wrap_model_call(
 self,
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
 ) -> ModelResponse:
 for attempt in range(self.max_retries):
 try:
 return handler(request)
 except Exception as e:
 if attempt == self.max_retries - 1:
 raise
 print(f"Retry {attempt + 1}/{self.max_retries} after error: {e}")

```

## [​](#state-updates)State updates

Both node-style and wrap-style hooks can update agent state. The mechanism differs:

{' ' * (self.list_depth - 1)}- **Node-style hooks** (`before_agent`, `before_model`, `after_model`, `after_agent`): Return a dict directly. The dict is applied to the agent state using the graph’s reducers.

{' ' * (self.list_depth - 1)}- **Wrap-style hooks** (`wrap_model_call`, `wrap_tool_call`): For model calls, return [`ExtendedModelResponse`](https://reference.langchain.com/python/langchain/agents/middleware/types/ExtendedModelResponse) with a [`Command`](https://reference.langchain.com/python/langgraph/types/Command) to inject state updates alongside the model response. For tool calls, return a [`Command`](https://reference.langchain.com/python/langgraph/types/Command) directly. Use these when you need to track or update state based on logic that runs during the model or tool call, such as summarization trigger points, usage metadata, or custom fields calculated from the request or response.

### [​](#node-style-hooks-2)Node-style hooks

Return a dict from a node-style hook to merge updates into agent state. The dict keys map to state fields.

```
from langchain.agents.middleware import after_model, AgentState
from langgraph.runtime import Runtime
from typing import Any
from typing_extensions import NotRequired

class TrackingState(AgentState):
 model_call_count: NotRequired[int]

@after_model(state_schema=TrackingState)
def increment_after_model(state: TrackingState, runtime: Runtime) -> dict[str, Any] | None:
 return {"model_call_count": state.get("model_call_count", 0) + 1}

```

### [​](#wrap-style-hooks-2)Wrap-style hooks

Return a [`ExtendedModelResponse`](https://reference.langchain.com/python/langchain/agents/middleware/types/ExtendedModelResponse) with a [`Command`](https://reference.langchain.com/python/langgraph/types/Command) from `wrap_model_call` to inject state updates from the model call layer:

```
from typing import Callable
from langchain.agents.middleware import (
 wrap_model_call,
 ModelRequest,
 ModelResponse,
 AgentState,
 ExtendedModelResponse
)
from langgraph.types import Command
from typing_extensions import NotRequired

class UsageTrackingState(AgentState):
 """Agent state with token usage tracking."""

 last_model_call_tokens: NotRequired[int]

@wrap_model_call(state_schema=UsageTrackingState)
def track_usage(
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
) -> ExtendedModelResponse:
 response = handler(request)
 return ExtendedModelResponse(
 model_response=response,
 command=Command(update={"last_model_call_tokens": 150}),
 )

```

The [`Command`](https://reference.langchain.com/python/langgraph/types/Command) flows through the graph’s reducers, so updates are applied correctly and messages are additive rather than replacing existing state.

#### [​](#composition-with-multiple-middleware)Composition with multiple middleware

When multiple middleware layers return `ExtendedModelResponse`, their commands compose:

{' ' * (self.list_depth - 1)}- **Commands are applied through reducers:** Each `Command` becomes a separate state update. For messages, this means they are additive.

{' ' * (self.list_depth - 1)}- **Outer wins on conflicts:** For non-reducer state fields, commands are applied inner-first, then outer. The outermost middleware’s value takes precedence on conflicting keys.

{' ' * (self.list_depth - 1)}- **Retry-safe:** If the outer middleware implements logic that can result in multiple calls to `handler()` again (for example, retry logic), commands from earlier calls are discarded.

```
from typing import Annotated, Callable

from langchain.agents.middleware import (
 AgentMiddleware,
 AgentState,
 ExtendedModelResponse,
 ModelRequest,
 ModelResponse,
)
from langchain.messages import SystemMessage
from langgraph.types import Command
from typing_extensions import NotRequired

def _last_wins(_a: str, b: str) -> str:
 """Reducer: last writer wins (outer overwrites inner)."""
 return b

class CustomMiddlewareState(AgentState):
 """Agent state: trace_layer uses last-wins (outer wins), messages use additive reducer."""

 # Non-reducer field with last-wins: both middleware write; outermost value wins
 trace_layer: NotRequired[Annotated[str, _last_wins]]

class OuterMiddleware(AgentMiddleware):
 def wrap_model_call(
 self,
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
 ) -> ExtendedModelResponse:
 response = handler(request)
 return ExtendedModelResponse(
 model_response=response,
 command=Command(update={
 "trace_layer": "outer",
 "messages": [SystemMessage(content="[Outer ran]")],
 }),
 )

class InnerMiddleware(AgentMiddleware):
 """Adds trace_layer and message. Outer adds to same keys; trace_layer: outer wins, messages: additive."""

 def wrap_model_call(
 self,
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
 ):
 response = handler(request)
 return ExtendedModelResponse(
 model_response=response,
 command=Command(update={
 "trace_layer": "inner",
 "messages": [SystemMessage(content="[Inner ran]")],
 }),
 )

```

## [​](#create-middleware)Create middleware

You can create middleware in two ways:

## Decorator-based middleware
Quick and simple for single-hook middleware. Use decorators to wrap individual functions.

## Class-based middleware
More powerful for complex middleware with multiple hooks or configuration.

### [​](#decorator-based-middleware)Decorator-based middleware

Quick and simple for single-hook middleware. Use decorators to wrap individual functions.
**Available decorators:**
**Node-style:**

{' ' * (self.list_depth - 1)}- [`@before_agent`](https://reference.langchain.com/python/langchain/agents/middleware/types/before_agent) - Runs before agent starts (once per invocation)

{' ' * (self.list_depth - 1)}- [`@before_model`](https://reference.langchain.com/python/langchain/agents/middleware/types/before_model) - Runs before each model call

{' ' * (self.list_depth - 1)}- [`@after_model`](https://reference.langchain.com/python/langchain/agents/middleware/types/after_model) - Runs after each model response

{' ' * (self.list_depth - 1)}- [`@after_agent`](https://reference.langchain.com/python/langchain/agents/middleware/types/after_agent) - Runs after agent completes (once per invocation)

**Wrap-style:**

{' ' * (self.list_depth - 1)}- [`@wrap_model_call`](https://reference.langchain.com/python/langchain/agents/middleware/types/wrap_model_call) - Wraps each model call with custom logic

{' ' * (self.list_depth - 1)}- [`@wrap_tool_call`](https://reference.langchain.com/python/langchain/agents/middleware/types/wrap_tool_call) - Wraps each tool call with custom logic

**Convenience:**

{' ' * (self.list_depth - 1)}- [`@dynamic_prompt`](https://reference.langchain.com/python/langchain/agents/middleware/types/dynamic_prompt) - Generates dynamic system prompts

**Example:**

```
from langchain.agents.middleware import (
 before_model,
 wrap_model_call,
 AgentState,
 ModelRequest,
 ModelResponse,
)
from langchain.agents import create_agent
from langgraph.runtime import Runtime
from typing import Any, Callable

@before_model
def log_before_model(state: AgentState, runtime: Runtime) -> dict[str, Any] | None:
 print(f"About to call model with {len(state['messages'])} messages")
 return None

@wrap_model_call
def retry_model(
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
) -> ModelResponse:
 for attempt in range(3):
 try:
 return handler(request)
 except Exception as e:
 if attempt == 2:
 raise
 print(f"Retry {attempt + 1}/3 after error: {e}")

agent = create_agent(
 model="gpt-5.4",
 middleware=[log_before_model, retry_model],
 tools=[...],
)

```

**When to use decorators:**

{' ' * (self.list_depth - 1)}- Single hook needed

{' ' * (self.list_depth - 1)}- No complex configuration

{' ' * (self.list_depth - 1)}- Quick prototyping

### [​](#class-based-middleware)Class-based middleware

More powerful for complex middleware with multiple hooks or configuration. Use classes when you need to define both sync and async implementations for the same hook, or when you want to combine multiple hooks in a single middleware.
**Example:**

```
from langchain.agents.middleware import (
 AgentMiddleware,
 AgentState,
 ModelRequest,
 ModelResponse,
)
from langgraph.runtime import Runtime
from typing import Any, Callable

class LoggingMiddleware(AgentMiddleware):
 def before_model(self, state: AgentState, runtime: Runtime) -> dict[str, Any] | None:
 print(f"About to call model with {len(state['messages'])} messages")
 return None

 def after_model(self, state: AgentState, runtime: Runtime) -> dict[str, Any] | None:
 print(f"Model returned: {state['messages'][-1].content}")
 return None

 async def abefore_model(
 self, state: AgentState, runtime: Runtime
 ) -> dict[str, Any] | None:
 # Async version of before_model
 return None

 async def aafter_model(
 self, state: AgentState, runtime: Runtime
 ) -> dict[str, Any] | None:
 # Async version of after_model
 print(f"Model returned: {state['messages'][-1].content}")
 return None

agent = create_agent(
 model="gpt-5.4",
 middleware=[LoggingMiddleware()],
 tools=[...],
)

```

**When to use classes:**

{' ' * (self.list_depth - 1)}- Defining both sync and async implementations for the same hook

{' ' * (self.list_depth - 1)}- Multiple hooks needed in a single middleware

{' ' * (self.list_depth - 1)}- Complex configuration required (e.g., configurable thresholds, custom models)

{' ' * (self.list_depth - 1)}- Reuse across projects with init-time configuration

## [​](#custom-state-schema)Custom state schema

If your middleware needs to track state across hooks, middleware can extend the agent’s state with custom properties. This enables middleware to:

{' ' * (self.list_depth - 1)}- 
**Track state across execution**: Maintain counters, flags, or other values that persist throughout the agent’s execution lifecycle

{' ' * (self.list_depth - 1)}- 
**Share data between hooks**: Pass information from `before_model` to `after_model` or between different middleware instances

{' ' * (self.list_depth - 1)}- 
**Implement cross-cutting concerns**: Add functionality like rate limiting, usage tracking, user context, or audit logging without modifying the core agent logic

{' ' * (self.list_depth - 1)}- 
**Make conditional decisions**: Use accumulated state to determine whether to continue execution, jump to different nodes, or modify behavior dynamically

{' ' * (self.list_depth - 1)}- Decorator
{' ' * (self.list_depth - 1)}- Class
```
from langchain.agents import create_agent
from langchain.messages import HumanMessage
from langchain.agents.middleware import AgentState, before_model, after_model
from typing_extensions import NotRequired
from typing import Any
from langgraph.runtime import Runtime

class CustomState(AgentState):
 model_call_count: NotRequired[int]
 user_id: NotRequired[str]

@before_model(state_schema=CustomState, can_jump_to=["end"])
def check_call_limit(state: CustomState, runtime: Runtime) -> dict[str, Any] | None:
 count = state.get("model_call_count", 0)
 if count > 10:
 return {"jump_to": "end"}
 return None

@after_model(state_schema=CustomState)
def increment_counter(state: CustomState, runtime: Runtime) -> dict[str, Any] | None:
 return {"model_call_count": state.get("model_call_count", 0) + 1}

agent = create_agent(
 model="gpt-5.4",
 middleware=[check_call_limit, increment_counter],
 tools=[],
)

# Invoke with custom state
result = agent.invoke({
 "messages": [HumanMessage("Hello")],
 "model_call_count": 0,
 "user_id": "user-123",
})

```

```
from langchain.agents import create_agent
from langchain.messages import HumanMessage
from langchain.agents.middleware import AgentState, AgentMiddleware
from typing_extensions import NotRequired
from typing import Any

class CustomState(AgentState):
 model_call_count: NotRequired[int]
 user_id: NotRequired[str]

class CallCounterMiddleware(AgentMiddleware[CustomState]):
 state_schema = CustomState

 def before_model(self, state: CustomState, runtime) -> dict[str, Any] | None:
 count = state.get("model_call_count", 0)
 if count > 10:
 return {"jump_to": "end"}
 return None

 def after_model(self, state: CustomState, runtime) -> dict[str, Any] | None:
 return {"model_call_count": state.get("model_call_count", 0) + 1}

agent = create_agent(
 model="gpt-5.4",
 middleware=[CallCounterMiddleware()],
 tools=[],
)

# Invoke with custom state
result = agent.invoke({
 "messages": [HumanMessage("Hello")],
 "model_call_count": 0,
 "user_id": "user-123",
})

```

## [​](#execution-order)Execution order

When using multiple middleware, understand how they execute:

```
agent = create_agent(
 model="gpt-5.4",
 middleware=[middleware1, middleware2, middleware3],
 tools=[...],
)

```

Execution flow**Before hooks run in order:**

{' ' * (self.list_depth - 1)}- `middleware1.before_agent()`

{' ' * (self.list_depth - 1)}- `middleware2.before_agent()`

{' ' * (self.list_depth - 1)}- `middleware3.before_agent()`
**Agent loop starts**

{' ' * (self.list_depth - 1)}- `middleware1.before_model()`

{' ' * (self.list_depth - 1)}- `middleware2.before_model()`

{' ' * (self.list_depth - 1)}- `middleware3.before_model()`
**Wrap hooks nest like function calls:**

{' ' * (self.list_depth - 1)}- `middleware1.wrap_model_call()` → `middleware2.wrap_model_call()` → `middleware3.wrap_model_call()` → model
**After hooks run in reverse order:**

{' ' * (self.list_depth - 1)}- `middleware3.after_model()`

{' ' * (self.list_depth - 1)}- `middleware2.after_model()`

{' ' * (self.list_depth - 1)}- `middleware1.after_model()`
**Agent loop ends**

{' ' * (self.list_depth - 1)}- `middleware3.after_agent()`

{' ' * (self.list_depth - 1)}- `middleware2.after_agent()`

{' ' * (self.list_depth - 1)}- `middleware1.after_agent()`

**Key rules:**

{' ' * (self.list_depth - 1)}- `before_*` hooks: First to last

{' ' * (self.list_depth - 1)}- `after_*` hooks: Last to first (reverse)

{' ' * (self.list_depth - 1)}- `wrap_*` hooks: Nested (first middleware wraps all others)

## [​](#agent-jumps)Agent jumps

To exit early from middleware, return a dictionary with `jump_to`:
**Available jump targets:**

{' ' * (self.list_depth - 1)}- `'end'`: Jump to the end of the agent execution (or the first `after_agent` hook)

{' ' * (self.list_depth - 1)}- `'tools'`: Jump to the tools node

{' ' * (self.list_depth - 1)}- `'model'`: Jump to the model node (or the first `before_model` hook)

{' ' * (self.list_depth - 1)}- Decorator
{' ' * (self.list_depth - 1)}- Class
```
from langchain.agents.middleware import after_model, hook_config, AgentState
from langchain.messages import AIMessage
from langgraph.runtime import Runtime
from typing import Any

@after_model
@hook_config(can_jump_to=["end"])
def check_for_blocked(state: AgentState, runtime: Runtime) -> dict[str, Any] | None:
 last_message = state["messages"][-1]
 if "BLOCKED" in last_message.content:
 return {
 "messages": [AIMessage("I cannot respond to that request.")],
 "jump_to": "end"
 }
 return None

```

```
from langchain.agents.middleware import AgentMiddleware, hook_config, AgentState
from langchain.messages import AIMessage
from langgraph.runtime import Runtime
from typing import Any

class BlockedContentMiddleware(AgentMiddleware):
 @hook_config(can_jump_to=["end"])
 def after_model(self, state: AgentState, runtime: Runtime) -> dict[str, Any] | None:
 last_message = state["messages"][-1]
 if "BLOCKED" in last_message.content:
 return {
 "messages": [AIMessage("I cannot respond to that request.")],
 "jump_to": "end"
 }
 return None

```

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- Keep middleware focused - each should do one thing well

{' ' * (self.list_depth - 1)}- Handle errors gracefully - don’t let middleware errors crash the agent

{' ' * (self.list_depth - 1)}- **Use appropriate hook types**:

{' ' * (self.list_depth - 1)}- Node-style for sequential logic (logging, validation)

{' ' * (self.list_depth - 1)}- Wrap-style for control flow (retry, fallback, caching)

{' ' * (self.list_depth - 1)}- Clearly document any custom state properties

{' ' * (self.list_depth - 1)}- Unit test middleware independently before integrating

{' ' * (self.list_depth - 1)}- Consider execution order - place critical middleware first in the list

{' ' * (self.list_depth - 1)}- Use built-in middleware when possible

## [​](#examples)Examples

### [​](#dynamic-prompt)Dynamic prompt

Dynamically modify the system prompt at runtime to inject context, user-specific instructions, or other information before each model call. This is one of the most common middleware use cases.
Use the `system_message` field on `ModelRequest` to read and modify the system prompt. It contains a [`SystemMessage`](https://reference.langchain.com/python/langchain-core/messages/system/SystemMessage) object (even if the agent was created with a string `system_prompt`).

{' ' * (self.list_depth - 1)}- Decorator
{' ' * (self.list_depth - 1)}- Class
```
from collections.abc import Callable

from langchain.agents.middleware import ModelRequest, ModelResponse, wrap_model_call
from langchain.messages import SystemMessage

@wrap_model_call
def add_context(
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
) -> ModelResponse:
 new_content = list(request.system_message.content_blocks) + [
 {"type": "text", "text": "Additional context."}
 ]
 new_system_message = SystemMessage(content=new_content)
 return handler(request.override(system_message=new_system_message))

```

```
from collections.abc import Callable

from langchain.agents.middleware import AgentMiddleware, ModelRequest, ModelResponse

class ContextMiddleware(AgentMiddleware):
 def wrap_model_call(
 self,
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
 ) -> ModelResponse:
 new_content = list(request.system_message.content_blocks) + [
 {"type": "text", "text": "Additional context."}
 ]
 new_system_message = SystemMessage(content=new_content)
 return handler(request.override(system_message=new_system_message))

```

{' ' * (self.list_depth - 1)}- `ModelRequest.system_message` is always a [`SystemMessage`](https://reference.langchain.com/python/langchain-core/messages/system/SystemMessage) object, even if the agent was created with `system_prompt="string"`

{' ' * (self.list_depth - 1)}- Use `SystemMessage.content_blocks` to access content as a list of blocks, regardless of whether the original content was a string or list

{' ' * (self.list_depth - 1)}- When modifying system messages, use `content_blocks` and append new blocks to preserve existing structure

{' ' * (self.list_depth - 1)}- You can pass [`SystemMessage`](https://reference.langchain.com/python/langchain-core/messages/system/SystemMessage) objects directly to `create_agent`’s `system_prompt` parameter for advanced use cases like cache control

### [​](#dynamic-model-selection)Dynamic model selection

{' ' * (self.list_depth - 1)}- Decorator
{' ' * (self.list_depth - 1)}- Class
```
from collections.abc import Callable

from langchain.agents.middleware import ModelRequest, ModelResponse, wrap_model_call
from langchain.chat_models import init_chat_model

complex_model = init_chat_model("claude-sonnet-4-6")
simple_model = init_chat_model("claude-haiku-4-5-20251001")

@wrap_model_call
def dynamic_model(
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
) -> ModelResponse:
 if len(request.messages) > 10:
 model = complex_model
 else:
 model = simple_model
 return handler(request.override(model=model))

```

```
from collections.abc import Callable

from langchain.agents.middleware import AgentMiddleware, ModelRequest, ModelResponse
from langchain.chat_models import init_chat_model

complex_model = init_chat_model("claude-sonnet-4-6")
simple_model = init_chat_model("claude-haiku-4-5-20251001")

class DynamicModelMiddleware(AgentMiddleware):
 def wrap_model_call(
 self,
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
 ) -> ModelResponse:
 if len(request.messages) > 10:
 model = complex_model
 else:
 model = simple_model
 return handler(request.override(model=model))

```

### [​](#dynamically-selecting-tools)Dynamically selecting tools

Select relevant tools at runtime to improve performance and accuracy. This section covers filtering pre-registered tools. For registering tools that are discovered at runtime (e.g., from MCP servers), see [Runtime tool registration](/oss/python/langchain/tools#dynamic-tool-selection).
**Benefits:**

{' ' * (self.list_depth - 1)}- **Shorter prompts** - Reduce complexity by exposing only relevant tools

{' ' * (self.list_depth - 1)}- **Better accuracy** - Models choose correctly from fewer options

{' ' * (self.list_depth - 1)}- **Permission control** - Dynamically filter tools based on user access

{' ' * (self.list_depth - 1)}- Decorator
{' ' * (self.list_depth - 1)}- Class
```
from langchain.agents import create_agent
from langchain.agents.middleware import wrap_model_call, ModelRequest, ModelResponse
from typing import Callable

@wrap_model_call
def select_tools(
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
) -> ModelResponse:
 """Middleware to select relevant tools based on state/context."""
 # Select a small, relevant subset of tools based on state/context
 relevant_tools = select_relevant_tools(request.state, request.runtime)
 return handler(request.override(tools=relevant_tools))

agent = create_agent(
 model="gpt-5.4",
 tools=all_tools, # All available tools need to be registered upfront
 middleware=[select_tools],
)

```

```
from langchain.agents import create_agent
from langchain.agents.middleware import AgentMiddleware, ModelRequest, ModelResponse
from typing import Callable

class ToolSelectorMiddleware(AgentMiddleware):
 def wrap_model_call(
 self,
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
 ) -> ModelResponse:
 """Middleware to select relevant tools based on state/context."""
 # Select a small, relevant subset of tools based on state/context
 relevant_tools = select_relevant_tools(request.state, request.runtime)
 return handler(request.override(tools=relevant_tools))

agent = create_agent(
 model="gpt-5.4",
 tools=all_tools, # All available tools need to be registered upfront
 middleware=[ToolSelectorMiddleware()],
)

```

### [​](#tool-call-monitoring)Tool call monitoring

{' ' * (self.list_depth - 1)}- Decorator
{' ' * (self.list_depth - 1)}- Class
```
from collections.abc import Callable

from langchain.agents.middleware import wrap_tool_call
from langchain.messages import ToolMessage
from langchain.tools.tool_node import ToolCallRequest
from langgraph.types import Command

@wrap_tool_call
def monitor_tool(
 request: ToolCallRequest,
 handler: Callable[[ToolCallRequest], ToolMessage | Command],
) -> ToolMessage | Command:
 print(f"Executing tool: {request.tool_call['name']}")
 print(f"Arguments: {request.tool_call['args']}")
 try:
 result = handler(request)
 print("Tool completed successfully")
 return result
 except Exception as e:
 print(f"Tool failed: {e}")
 raise

```

```
from collections.abc import Callable

from langchain.agents.middleware import AgentMiddleware
from langchain.messages import ToolMessage
from langchain.tools.tool_node import ToolCallRequest
from langgraph.types import Command

class ToolMonitoringMiddleware(AgentMiddleware):
 def wrap_tool_call(
 self,
 request: ToolCallRequest,
 handler: Callable[[ToolCallRequest], ToolMessage | Command],
 ) -> ToolMessage | Command:
 print(f"Executing tool: {request.tool_call['name']}")
 print(f"Arguments: {request.tool_call['args']}")
 try:
 result = handler(request)
 print("Tool completed successfully")
 return result
 except Exception as e:
 print(f"Tool failed: {e}")
 raise

```

### [​](#prompt-caching-anthropic)Prompt caching (Anthropic)

When working with Anthropic models, use structured content blocks with cache control directives to cache large system prompts:

{' ' * (self.list_depth - 1)}- Decorator
{' ' * (self.list_depth - 1)}- Class
```
from langchain.agents.middleware import wrap_model_call, ModelRequest, ModelResponse
from langchain.messages import SystemMessage
from typing import Callable

@wrap_model_call
def add_cached_context(
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
) -> ModelResponse:
 # Always work with content blocks
 new_content = list(request.system_message.content_blocks) + [
 {
 "type": "text",
 "text": "Here is a large document to analyze:\n\n<document>...</document>",
 # content up until this point is cached
 "cache_control": {"type": "ephemeral"}
 }
 ]

 new_system_message = SystemMessage(content=new_content)
 return handler(request.override(system_message=new_system_message))

```

```
from langchain.agents.middleware import AgentMiddleware, ModelRequest, ModelResponse
from langchain.messages import SystemMessage
from typing import Callable

class CachedContextMiddleware(AgentMiddleware):
 def wrap_model_call(
 self,
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
 ) -> ModelResponse:
 # Always work with content blocks
 new_content = list(request.system_message.content_blocks) + [
 {
 "type": "text",
 "text": "Here is a large document to analyze:\n\n<document>...</document>",
 "cache_control": {"type": "ephemeral"} # This content will be cached
 }
 ]

 new_system_message = SystemMessage(content=new_content)
 return handler(request.override(system_message=new_system_message))

```

**Notes:**

{' ' * (self.list_depth - 1)}- `ModelRequest.system_message` is always a [`SystemMessage`](https://reference.langchain.com/python/langchain-core/messages/system/SystemMessage) object, even if the agent was created with `system_prompt="string"`

{' ' * (self.list_depth - 1)}- Use `SystemMessage.content_blocks` to access content as a list of blocks, regardless of whether the original content was a string or list

{' ' * (self.list_depth - 1)}- When modifying system messages, use `content_blocks` and append new blocks to preserve existing structure

{' ' * (self.list_depth - 1)}- You can pass [`SystemMessage`](https://reference.langchain.com/python/langchain-core/messages/system/SystemMessage) objects directly to `create_agent`’s `system_prompt` parameter for advanced use cases like cache control

:::

## [​](#additional-resources)Additional resources

{' ' * (self.list_depth - 1)}- [Middleware API reference](https://reference.langchain.com/python/langchain/middleware/)

{' ' * (self.list_depth - 1)}- [Built-in middleware](/oss/python/langchain/middleware/built-in)

{' ' * (self.list_depth - 1)}- [Testing agents](/oss/python/langchain/test)

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/middleware/custom.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Prebuilt middlewarePrevious](/oss/python/langchain/middleware/built-in)[OverviewNext](/oss/python/langchain/frontend/overview)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
