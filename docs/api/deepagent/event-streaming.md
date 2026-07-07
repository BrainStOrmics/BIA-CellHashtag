# Event streaming

[Core capabilities](/oss/python/deepagents/harness)

# Event streaming
Copy page

Stream subagents, messages, tool calls, and final output from Deep Agents.Copy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.This page covers streaming concerns specific to Deep Agents—most importantly, streaming from delegated subagents via `stream.subagents`. For general agent streaming (`stream.messages`, `stream.values`, tool calls, custom updates), see [LangChain Event Streaming](/oss/python/langchain/event-streaming).

## [​](#stream-subagents)Stream subagents

Deep Agents add a subagent projection on top of LangGraph streaming. Use `stream.subagents` when you want one stream handle per delegated `task` call. The projection is lightweight: it discovers subagent tasks first, and message, tool-call, and value streams are opened only when you access them on a subagent handle.

```
stream = agent.stream_events({
 "messages": [{"role": "user", "content": "Write me a haiku about the sea"}],
}, version="v3")

for subagent in stream.subagents:
 print(subagent.name, subagent.path, subagent.status)

 for message in subagent.messages:
 print(message.text)

```

## [​](#subagent-stream-fields)Subagent stream fields

Each subagent stream exposes the same kinds of projections as the parent run, such as messages, tool calls, nested subagents, and final output. For the general parent-run streaming model, see [LangChain Event Streaming](/oss/python/langchain/event-streaming).
Python uses snake_case projection names such as `tool_calls`. Each subagent stream can expose `.messages`, `.tool_calls`, `.values`, `.subagents`, and `.output`.

| 
 | Field | Description
 | `name` | Subagent name.
 | `messages` | Messages emitted by the subagent.
 | `subagents` | Nested subagent invocations.
 | `output` | Final subagent state, or completion signal for the delegated task.
 | `path` | Namespace path for the subagent stream.
 | `status` | Lifecycle status such as `started`, `completed`, `failed`, or `interrupted`.
 | `tool_calls` | Tool calls scoped to the subagent.

## [​](#track-subagent-lifecycle)Track subagent lifecycle

Use `stream.subagents` when you only need to show which subagents started and finished. You do not need to subscribe to message or value streams unless you access those projections on an individual subagent.

```
stream = agent.stream_events(input, version="v3")

running = 0
completed = 0
failed = 0

for subagent in stream.subagents:
 running += 1
 print(f"{subagent.name}: started")

 try:
 _ = subagent.output
 running -= 1
 completed += 1
 print(f"{subagent.name}: completed")
 except Exception:
 running -= 1
 failed += 1
 print(f"{subagent.name}: failed")

```

## [​](#stream-messages)Stream messages

Deep Agents can emit messages from the coordinator agent and from delegated subagents. Use `stream.messages` for top-level messages and `subagent.messages` for each delegated subagent.

```
stream = agent.stream_events(input, version="v3")

for message in stream.messages:
 print("[coordinator]", message.text)

for subagent in stream.subagents:
 for message in subagent.messages:
 print(f"[{subagent.name}]", message.text)

```

## [​](#stream-tool-calls)Stream tool calls

Deep Agents expose tool calls at each level of the agent tree. Use the top-level `stream.tool_calls` for coordinator tools and each `subagent.tool_calls` for delegated work.

```
stream = agent.stream_events(input, version="v3")

for call in stream.tool_calls:
 print("[coordinator tool]", call.tool_name, call.input)
 print(call.completed, call.error)

for subagent in stream.subagents:
 for call in subagent.tool_calls:
 print(f"[{subagent.name} tool]", call.tool_name, call.input)
 for delta in call.output_deltas:
 print(delta, end="", flush=True)

 if call.completed and call.error is None:
 print(call.output)
 elif call.error is not None:
 print(call.error)

```

## [​](#stream-nested-work)Stream nested work

You can recurse into a subagent stream to observe nested subagents, messages, and tool calls.

```
stream = agent.stream_events(input, version="v3")

for subagent in stream.subagents:
 print(f"subagent {subagent.name}: {subagent.status}")

 for tool_call in subagent.tool_calls:
 print(f"{tool_call.tool_name}({tool_call.input})")
 for delta in tool_call.output_deltas:
 print(delta, end="", flush=True)

 for nested in subagent.subagents:
 print(f"nested subagent {nested.name}: {nested.status}")

```

## [​](#consume-concurrently)Consume concurrently

Coordinator and subagent output often interleave. Consume projections concurrently when you need live UI updates.
For concurrent consumption in async code, use `astream_events` with `asyncio.gather`:

```
import asyncio

stream = await agent.astream_events(input, version="v3")

async def consume_coordinator():
 async for message in stream.messages:
 print("[coordinator]", await message.text)

async def consume_subagents():
 async for subagent in stream.subagents:
 async for message in subagent.messages:
 print(f"[{subagent.name}]", await message.text)

await asyncio.gather(consume_coordinator(), consume_subagents())

```

For synchronous code, use `stream.interleave(...)` instead:

```
stream = agent.stream_events(input, version="v3")

for name, item in stream.interleave("messages", "subagents"):
 if name == "messages":
 print("[coordinator]", item.text)
 else:
 for message in item.messages:
 print(f"[{item.name}]", message.text)

```

When you need exact arrival order across the coordinator and all subagents, iterate raw protocol events and use `namespace` to identify the source:

```
stream = agent.stream_events(input, version="v3")

for event in stream:
 if event.get("method") != "messages":
 continue

 payload = event["params"]["data"][0]
 if not isinstance(payload, dict):
 continue
 if payload.get("event") != "content-block-delta":
 continue

 block = payload.get("delta") or {}
 if block.get("type") == "text-delta":
 source = "subagent" if event["params"]["namespace"] else "coordinator"
 print(f"[{source}] {block['text']}")

```

## [​](#subagents-versus-subgraphs)Subagents versus subgraphs

`stream.subgraphs` shows graph execution structure. `stream.subagents` shows product-level Deep Agents task delegations. Use `stream.subagents` for user-facing UI because it hides internal graph nodes and exposes the subagent concept directly.

## [​](#related)Related

{' ' * (self.list_depth - 1)}- [LangChain Event Streaming](/oss/python/langchain/event-streaming) covers general agent message and tool-call streaming concepts.

{' ' * (self.list_depth - 1)}- [Subagent frontend streaming](/oss/python/deepagents/frontend/subagent-streaming) shows UI patterns that separate coordinator messages from subagent cards.

{' ' * (self.list_depth - 1)}- [LangGraph Event Streaming](/oss/python/langgraph/event-streaming) covers the underlying graph streaming model.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/event-streaming.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[ProfilesPrevious](/oss/python/deepagents/profiles)[StreamingNext](/oss/python/deepagents/streaming)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
