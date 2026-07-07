# Event streaming

[Core components](/oss/python/langchain/agents)

# Event streaming
Copy page

Stream real-time updates from LangChain agent runsCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.LangChain agents are built on LangGraph, so they support the same streaming stack with agent-focused projections for messages, tool calls, state, and custom updates.
For most application and frontend use cases, use **Event Streaming** through `stream_events(..., version="v3")`. Event Streaming returns a run object with typed projections, so each projection can be consumed independently instead of parsing stream-mode tuples.

```
from langchain.agents import create_agent

def get_weather(city: str) -> str:
 """Get weather for a city."""
 return f"It's always sunny in {city}!"

agent = create_agent(
 model="gpt-5-nano",
 tools=[get_weather],
)

stream = agent.stream_events({
 "messages": [{"role": "user", "content": "What is the weather in SF?"}],
}, version="v3")

for message in stream.messages:
 for delta in message.text:
 print(delta, end="", flush=True)

final_state = stream.output

```

## [​](#what-you-can-stream)What you can stream

| 
 | Projection | Use
 | `for event in stream` | Raw protocol events with full envelope and access to every channel.
 | `stream.messages` | Model message streams, one per LLM call.
 | `message.text` | Text deltas and final text for a message.
 | `message.reasoning` | Reasoning deltas for models that expose reasoning content.
 | `message.tool_calls` | Tool-call argument chunks and finalized tool calls.
 | `message.output` | Final message object after the model call completes.
 | `stream.values` | Agent state snapshots.
 | `stream.output` | Final agent state.
 | `stream.subgraphs` | Nested graph runs (sub-agents and plain subgraphs).
 | `stream.extensions` | Custom transformer projections.
 | `stream.tool_calls` | Tool execution lifecycle, inputs, output deltas, final output, and errors.
`stream.messages` yields `ChatModelStream` objects. Each message stream exposes `.text`, `.reasoning`, `.tool_calls`, and `.output`. Sync projections are iterable for live deltas and drainable for final values: use `str(message.text)` for final text and `message.tool_calls.get()` for finalized tool calls.

## [​](#agent-messages)Agent messages

Use `stream.messages` when you want model output from each LLM call.

```
stream = agent.stream_events(input, version="v3")

for message in stream.messages:
 print(f"[{message.node}] ", end="")
 for delta in message.text:
 print(delta, end="", flush=True)

 full_message = message.output
 usage = full_message.usage_metadata
 if usage:
 print(usage)

```

`message.output` gives you the finalized AI message, including provider-specific content blocks. In TypeScript, use `message.usage` when you only need token counts or other usage metadata; in Python, read usage from `message.output.usage_metadata`.

## [​](#reasoning-content)Reasoning content

Reasoning content uses the same shape as text content, but it is available only when the selected model emits reasoning blocks.

```
stream = agent.stream_events(input, version="v3")

for message in stream.messages:
 for delta in message.reasoning:
 print(f"[thinking] {delta}", end="", flush=True)

 for delta in message.text:
 print(delta, end="", flush=True)

```

See the [reasoning guide](/oss/python/langchain/models#reasoning) and your provider’s integration page for model configuration details.

## [​](#tool-calls)Tool calls

There are two useful tool-call projections:

{' ' * (self.list_depth - 1)}- `message.tool_calls` streams tool-call argument chunks while the model is producing the tool call.

{' ' * (self.list_depth - 1)}- `stream.tool_calls` streams the lifecycle of tool execution after the tool call starts.

```
stream = agent.stream_events(input, version="v3")

for message in stream.messages:
 for chunk in message.tool_calls:
 print(f"tool call chunk: {chunk}")

 finalized = message.tool_calls.get()
 if finalized:
 print(f"finalized tool calls: {finalized}")

for call in stream.tool_calls:
 print(f"{call.tool_name}({call.input})")
 for delta in call.output_deltas:
 print(delta, end="", flush=True)
 print(call.output, call.error)

```

## [​](#streaming-sub-agents)Streaming sub-agents

When a `create_agent` call invokes another `create_agent` (via a wrapping tool, typically), the inner agent’s events flow at a nested namespace and surface as a handle on `stream.subgraphs`. Each handle exposes the inner agent’s own `.messages`, `.values`, `.tool_calls`, and `.output` projections. The `name=` you pass to `create_agent` becomes `subagent.graph_name` (Python) / `subagent.name` (JS), which lets you filter and label per agent.
Every nested `CompiledStateGraph` shows up on `stream.subgraphs` — `create_agent` instances are one specific kind. Filter on the name to act only on the ones you care about.

```
from langchain.agents import create_agent
from langchain.chat_models import init_chat_model

def get_weather(city: str) -> str:
 """Get weather for a given city."""
 return f"It's always sunny in {city}!"

weather_agent = create_agent(
 model=init_chat_model("openai:gpt-5.4"),
 tools=[get_weather],
 name="weather_agent",
)

def call_weather(query: str) -> str:
 """Query the weather agent."""
 result = weather_agent.invoke({"messages": [{"role": "user", "content": query}]})
 return result["messages"][-1].text

supervisor = create_agent(
 model=init_chat_model("openai:gpt-5.4"),
 tools=[call_weather],
 name="supervisor",
)

stream = supervisor.stream_events(
 {"messages": [{"role": "user", "content": "What's the weather in Boston?"}]},
 version="v3",
)

for subagent in stream.subgraphs:
 if subagent.graph_name != "weather_agent":
 continue
 print(f"{subagent.graph_name}: ", end="")
 for message in subagent.messages:
 for token in message.text:
 print(token, end="", flush=True)
 print()

```

The same projection covers plain `StateGraph` subgraphs invoked from a tool — set `name=` on `.compile(name=...)` to get a label in `subagent.graph_name`. There’s no separate sub-agent-only projection; the filter is what you write into your loop.

## [​](#state-and-final-output)State and final output

Use `stream.values` for state snapshots and `stream.output` for the final agent state.

```
stream = agent.stream_events(input, version="v3")

for snapshot in stream.values:
 print(snapshot)

final_state = stream.output

```

## [​](#multiple-projections)Multiple projections

For concurrent consumption in async code, use `astream_events` with `asyncio.gather`:

```
import asyncio

stream = await agent.astream_events(input, version="v3")

async def consume_messages():
 async for message in stream.messages:
 print(await message.text)

async def consume_tool_calls():
 async for call in stream.tool_calls:
 print(call.tool_name, call.input)

await asyncio.gather(consume_messages(), consume_tool_calls())

```

For synchronous code, use `stream.interleave(...)` instead:

```
stream = agent.stream_events(input, version="v3")

for name, item in stream.interleave("messages", "tool_calls", "values"):
 if name == "messages":
 print(item.text)
 elif name == "tool_calls":
 print(item.tool_name, item.input)
 elif name == "values":
 print(item)

```

To access channels that aren’t exposed as typed projections, or to inspect the full event envelope, iterate raw protocol events:

```
for event in stream:
 print(event["method"], event["params"]["namespace"], event["params"]["data"])

```

## [​](#custom-updates)Custom updates

Use custom stream transformers when your application needs a projection that is not built in, such as retrieval progress, artifacts, or domain-specific events.

```
stream = agent.stream_events(
 input,
 version="v3",
 transformers=[ToolActivityTransformer],
)

for activity in stream.extensions["tool_activity"]:
 print(activity)

```

See [Build your own projection](/oss/python/langgraph/event-streaming#build-your-own-projection) for the transformer contract.

## [​](#related)Related

{' ' * (self.list_depth - 1)}- [Streaming](/oss/python/langchain/streaming) covers low-level Pregel stream modes.

{' ' * (self.list_depth - 1)}- [Build your own projection](/oss/python/langgraph/event-streaming#build-your-own-projection) covers writing application-specific projections.

{' ' * (self.list_depth - 1)}- [Frontend streaming patterns](/oss/python/langchain/frontend/overview) shows UI use cases built on streamed state.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/event-streaming.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Short-term memoryPrevious](/oss/python/langchain/short-term-memory)[StreamingNext](/oss/python/langchain/streaming)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
