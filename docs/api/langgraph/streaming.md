# Streaming

[Capabilities](/oss/python/langgraph/persistence)

# Streaming
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.For new applications, we recommend [event streaming](/oss/python/langgraph/event-streaming)—the typed-projection API introduced in LangGraph v1.2. Event streaming gives you separate iterators per projection (messages, values, subgraphs, output) so you can consume them independently instead of branching on `stream_mode` chunks.
This page covers LangGraph’s stream-mode API. It exposes graph execution through stream modes such as `updates`, `values`, `messages`, `custom`, `checkpoints`, `tasks`, and `debug`. Use it when you need direct access to graph-runtime events or specific stream-mode output.

## [​](#get-started)Get started

### [​](#basic-usage)Basic usage

LangGraph graphs expose the [`stream`](https://reference.langchain.com/python/langgraph/pregel/#langgraph.pregel.Pregel.stream) (sync) and [`astream`](https://reference.langchain.com/python/langgraph/pregel/#langgraph.pregel.Pregel.astream) (async) methods to yield streamed outputs as iterators. Pass one or more [stream modes](#stream-modes) to control what data you receive.

```
for chunk in graph.stream(
 {"topic": "ice cream"},
 stream_mode=["updates", "custom"],
 version="v2",
):
 if chunk["type"] == "updates":
 for node_name, state in chunk["data"].items():
 print(f"Node {node_name} updated: {state}")
 elif chunk["type"] == "custom":
 print(f"Status: {chunk['data']['status']}")

```

Output
```
Status: thinking of a joke...
Node generate_joke updated: {'joke': 'Why did the ice cream go to school? To get a sundae education!'}

```

Full example
```
from typing import TypedDict
from langgraph.graph import StateGraph, START, END
from langgraph.config import get_stream_writer

class State(TypedDict):
 topic: str
 joke: str

def generate_joke(state: State):
 writer = get_stream_writer()
 writer({"status": "thinking of a joke..."})
 return {"joke": f"Why did the {state['topic']} go to school? To get a sundae education!"}

graph = (
 StateGraph(State)
 .add_node(generate_joke)
 .add_edge(START, "generate_joke")
 .add_edge("generate_joke", END)
 .compile()
)

for chunk in graph.stream(
 {"topic": "ice cream"},
 stream_mode=["updates", "custom"],
 version="v2",
):
 if chunk["type"] == "updates":
 for node_name, state in chunk["data"].items():
 print(f"Node {node_name} updated: {state}")
 elif chunk["type"] == "custom":
 print(f"Status: {chunk['data']['status']}")

```
Output
```
Status: thinking of a joke...
Node generate_joke updated: {'joke': 'Why did the ice cream go to school? To get a sundae education!'}

```

Debug streaming events, inspect token-by-token LLM output, and monitor latency with [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langgraph-streaming). Follow the [tracing quickstart](/langsmith/trace-with-langgraph) to get set up.

### [​](#stream-output-format-v2)Stream output format (v2)

Requires LangGraph >= 1.1. All examples on this page use `version="v2"`.
Pass `version="v2"` to `stream()` or `astream()` to get a unified output format. Every chunk is a `StreamPart` dict with a consistent shape — regardless of stream mode, number of modes, or subgraph settings:

```
{
 "type": "values" | "updates" | "messages" | "custom" | "checkpoints" | "tasks" | "debug",
 "ns": (), # namespace tuple, populated for subgraph events
 "data": ..., # the actual payload (type varies by stream mode)
}

```

Each stream mode has a corresponding `TypedDict` containing [`ValuesStreamPart`](https://reference.langchain.com/python/langgraph/types/ValuesStreamPart), [`UpdatesStreamPart`](https://reference.langchain.com/python/langgraph/types/UpdatesStreamPart), [`MessagesStreamPart`](https://reference.langchain.com/python/langgraph/types/MessagesStreamPart), [`CustomStreamPart`](https://reference.langchain.com/python/langgraph/types/CustomStreamPart), [`CheckpointStreamPart`](https://reference.langchain.com/python/langgraph/types/CheckpointStreamPart), [`TasksStreamPart`](https://reference.langchain.com/python/langgraph/types/TasksStreamPart), [`DebugStreamPart`](https://reference.langchain.com/python/langgraph/types/DebugStreamPart). You can import these types from `langgraph.types`. The union type [`StreamPart`](https://reference.langchain.com/python/langgraph/types/StreamPart) is a disjoing union on `part["type"]`, enabling full type narrowing in editors and type checkers.
With v1 (default), the output format changes based on your streaming options (single mode returns raw data, multiple modes return `(mode, data)` tuples, subgraphs return `(namespace, data)` tuples). With v2, the format is always the same:
v2 (new)v1 (current default)
```
for chunk in graph.stream(inputs, stream_mode="updates", version="v2"):
 print(chunk["type"]) # "updates"
 print(chunk["ns"]) # ()
 print(chunk["data"]) # {"node_name": {"key": "value"}}

```

The v2 format also enables type narrowing, which means you can filter chunks by `chunk["type"]` and get the correct payload type. Each branch narrows `part["data"]` to the specific type for that mode:

```
for part in graph.stream(
 {"topic": "ice cream"},
 stream_mode=["values", "updates", "messages", "custom"],
 version="v2",
):
 if part["type"] == "values":
 # ValuesStreamPart — full state snapshot after each step
 print(f"State: topic={part['data']['topic']}")
 elif part["type"] == "updates":
 # UpdatesStreamPart — only the changed keys from each node
 for node_name, state in part["data"].items():
 print(f"Node `{node_name}` updated: {state}")
 elif part["type"] == "messages":
 # MessagesStreamPart — (message_chunk, metadata) from LLM calls
 msg, metadata = part["data"]
 print(msg.content, end="", flush=True)
 elif part["type"] == "custom":
 # CustomStreamPart — arbitrary data from get_stream_writer()
 print(f"Progress: {part['data']['progress']}%")

```

## [​](#stream-modes)Stream modes

Pass one or more of the following stream modes as a list to the [`stream`](https://reference.langchain.com/python/langgraph/graphs/#langgraph.graph.state.CompiledStateGraph.stream) or [`astream`](https://reference.langchain.com/python/langgraph/graphs/#langgraph.graph.state.CompiledStateGraph.astream) methods:

| 
 | Mode | Type | Description
 | [values](#graph-state) | [`ValuesStreamPart`](https://reference.langchain.com/python/langgraph/types/ValuesStreamPart) | Full state after each step.
 | [updates](#graph-state) | [`UpdatesStreamPart`](https://reference.langchain.com/python/langgraph/types/UpdatesStreamPart) | State updates after each step. Multiple updates in the same step are streamed separately.
 | [messages](#llm-tokens) | [`MessagesStreamPart`](https://reference.langchain.com/python/langgraph/types/MessagesStreamPart) | 2-tuples of (LLM token, metadata) from LLM calls.
 | [custom](#custom-data) | [`CustomStreamPart`](https://reference.langchain.com/python/langgraph/types/CustomStreamPart) | Custom data emitted from nodes via [`get_stream_writer`](https://reference.langchain.com/python/langgraph/config/get_stream_writer).
 | [checkpoints](#checkpoints) | [`CheckpointStreamPart`](https://reference.langchain.com/python/langgraph/types/CheckpointStreamPart) | Checkpoint events (same format as `get_state()`). Requires a checkpointer.
 | [tasks](#tasks) | [`TasksStreamPart`](https://reference.langchain.com/python/langgraph/types/TasksStreamPart) | Task start/finish events with results and errors. Requires a checkpointer.
 | [debug](#debug) | [`DebugStreamPart`](https://reference.langchain.com/python/langgraph/types/DebugStreamPart) | All available info — combines `checkpoints` and `tasks` with extra metadata.
[]()

### [​](#graph-state)Graph state

Use the stream modes `updates` and `values` to stream the state of the graph as it executes.

{' ' * (self.list_depth - 1)}- `updates` streams the **updates** to the state after each step of the graph.

{' ' * (self.list_depth - 1)}- `values` streams the **full value** of the state after each step of the graph.

```
from typing import TypedDict
from langgraph.graph import StateGraph, START, END

class State(TypedDict):
 topic: str
 joke: str

def refine_topic(state: State):
 return {"topic": state["topic"] + " and cats"}

def generate_joke(state: State):
 return {"joke": f"This is a joke about {state['topic']}"}

graph = (
 StateGraph(State)
 .add_node(refine_topic)
 .add_node(generate_joke)
 .add_edge(START, "refine_topic")
 .add_edge("refine_topic", "generate_joke")
 .add_edge("generate_joke", END)
 .compile()
)

```

{' ' * (self.list_depth - 1)}- updates
{' ' * (self.list_depth - 1)}- valuesUse this to stream only the **state updates** returned by the nodes after each step. The streamed outputs include the name of the node as well as the update.
```
for chunk in graph.stream(
 {"topic": "ice cream"},
 stream_mode="updates",
 version="v2",
):
 if chunk["type"] == "updates":
 for node_name, state in chunk["data"].items():
 print(f"Node `{node_name}` updated: {state}")

```
Output
```
Node `refine_topic` updated: {'topic': 'ice cream and cats'}
Node `generate_joke` updated: {'joke': 'This is a joke about ice cream and cats'}

```
Use this to stream the **full state** of the graph after each step.
```
for chunk in graph.stream(
 {"topic": "ice cream"},
 stream_mode="values",
 version="v2",
):
 if chunk["type"] == "values":
 print(f"topic: {chunk['data']['topic']}, joke: {chunk['data']['joke']}")

```
Output
```
topic: ice cream, joke:
topic: ice cream and cats, joke:
topic: ice cream and cats, joke: This is a joke about ice cream and cats

```

### [​](#llm-tokens)LLM tokens

Use the `messages` streaming mode to stream Large Language Model (LLM) outputs **token by token** from any part of your graph, including nodes, tools, subgraphs, or tasks.
The streamed output from [`messages` mode](#stream-modes) is a tuple `(message_chunk, metadata)` where:

{' ' * (self.list_depth - 1)}- `message_chunk`: the token or message segment from the LLM.

{' ' * (self.list_depth - 1)}- `metadata`: a dictionary containing details about the graph node and LLM invocation.

> 
If your LLM is not available as a LangChain integration, you can stream its outputs using `custom` mode instead. See [use with any LLM](#use-with-any-llm) for details.

**Manual config required for async in Python < 3.11**
When using Python < 3.11 with async code, you must explicitly pass [`RunnableConfig`](https://reference.langchain.com/python/langchain-core/runnables/config/RunnableConfig) to `ainvoke()` to enable proper streaming. See [Async with Python < 3.11](#async) for details or upgrade to Python 3.11+.

```
from dataclasses import dataclass

from langchain.chat_models import init_chat_model
from langgraph.graph import StateGraph, START

@dataclass
class MyState:
 topic: str
 joke: str = ""

model = init_chat_model(model="gpt-5.4-mini")

def call_model(state: MyState):
 """Call the LLM to generate a joke about a topic"""
 # Note that message events are emitted even when the LLM is run using .invoke rather than .stream
 model_response = model.invoke(
 [
 {"role": "user", "content": f"Generate a joke about {state.topic}"}
 ]
 )
 return {"joke": model_response.content}

graph = (
 StateGraph(MyState)
 .add_node(call_model)
 .add_edge(START, "call_model")
 .compile()
)

# The "messages" stream mode streams LLM tokens with metadata
# Use version="v2" for a unified StreamPart format
for chunk in graph.stream(
 {"topic": "ice cream"},
 stream_mode="messages",
 version="v2",
):
 if chunk["type"] == "messages":
 message_chunk, metadata = chunk["data"]
 if message_chunk.content:
 print(message_chunk.content, end="|", flush=True)

```

#### [​](#filter-by-llm-invocation)Filter by LLM invocation

You can associate `tags` with LLM invocations to filter the streamed tokens by LLM invocation.

```
from langchain.chat_models import init_chat_model

# model_1 is tagged with "joke"
model_1 = init_chat_model(model="gpt-5.4-mini", tags=['joke'])
# model_2 is tagged with "poem"
model_2 = init_chat_model(model="gpt-5.4-mini", tags=['poem'])

graph = ... # define a graph that uses these LLMs

# The stream_mode is set to "messages" to stream LLM tokens
# The metadata contains information about the LLM invocation, including the tags
async for chunk in graph.astream(
 {"topic": "cats"},
 stream_mode="messages",
 version="v2",
):
 if chunk["type"] == "messages":
 msg, metadata = chunk["data"]
 # Filter the streamed tokens by the tags field in the metadata to only include
 # the tokens from the LLM invocation with the "joke" tag
 if metadata["tags"] == ["joke"]:
 print(msg.content, end="|", flush=True)

```

Extended example: filtering by tags
```
from typing import TypedDict

from langchain.chat_models import init_chat_model
from langgraph.graph import START, StateGraph

# The joke_model is tagged with "joke"
joke_model = init_chat_model(model="gpt-5.4-mini", tags=["joke"])
# The poem_model is tagged with "poem"
poem_model = init_chat_model(model="gpt-5.4-mini", tags=["poem"])

class State(TypedDict):
 topic: str
 joke: str
 poem: str

async def call_model(state, config):
 topic = state["topic"]
 print("Writing joke...")
 # Note: Passing the config through explicitly is required for python < 3.11
 # Since context var support wasn't added before then: https://docs.python.org/3/library/asyncio-task.html#creating-tasks
 # The config is passed through explicitly to ensure the context vars are propagated correctly
 # This is required for Python < 3.11 when using async code. Please see the async section for more details
 joke_response = await joke_model.ainvoke(
 [{"role": "user", "content": f"Write a joke about {topic}"}],
 config,
 )
 print("\n\nWriting poem...")
 poem_response = await poem_model.ainvoke(
 [{"role": "user", "content": f"Write a short poem about {topic}"}],
 config,
 )
 return {"joke": joke_response.content, "poem": poem_response.content}

graph = (
 StateGraph(State)
 .add_node(call_model)
 .add_edge(START, "call_model")
 .compile()
)

# The stream_mode is set to "messages" to stream LLM tokens
# The metadata contains information about the LLM invocation, including the tags
async for chunk in graph.astream(
 {"topic": "cats"},
 stream_mode="messages",
 version="v2",
):
 if chunk["type"] == "messages":
 msg, metadata = chunk["data"]
 if metadata["tags"] == ["joke"]:
 print(msg.content, end="|", flush=True)

```

#### [​](#omit-messages-from-the-stream)Omit messages from the stream

Use the `nostream` tag to exclude LLM output from the stream entirely. Invocations tagged with `nostream` still run and produce output; their tokens are simply not emitted in `messages` mode.
This is useful when:

{' ' * (self.list_depth - 1)}- You need LLM output for internal processing (for example structured output) but do not want to stream it to the client

{' ' * (self.list_depth - 1)}- You stream the same content through a different channel (for example custom UI messages) and want to avoid duplicate output in the `messages` stream

```
from typing import Any, TypedDict

from langchain_anthropic import ChatAnthropic
from langgraph.graph import START, StateGraph

stream_model = ChatAnthropic(model_name="claude-haiku-4-5-20251001")
internal_model = ChatAnthropic(model_name="claude-haiku-4-5-20251001").with_config(
 {"tags": ["nostream"]}
)

class State(TypedDict):
 topic: str
 answer: str
 notes: str

def answer(state: State) -> dict[str, Any]:
 r = stream_model.invoke(
 [{"role": "user", "content": f"Reply briefly about {state['topic']}"}]
 )
 return {"answer": r.content}

def internal_notes(state: State) -> dict[str, Any]:
 # Tokens from this model are omitted from stream_mode="messages" because of nostream
 r = internal_model.invoke(
 [{"role": "user", "content": f"Private notes on {state['topic']}"}]
 )
 return {"notes": r.content}

graph = (
 StateGraph(State)
 .add_node("write_answer", answer)
 .add_node("internal_notes", internal_notes)
 .add_edge(START, "write_answer")
 .add_edge("write_answer", "internal_notes")
 .compile()
)

initial_state: State = {"topic": "AI", "answer": "", "notes": ""}
stream = graph.stream(initial_state, stream_mode="messages")

```

#### [​](#filter-by-node)Filter by node

To stream tokens only from specific nodes, use `stream_mode="messages"` and filter the outputs by the `langgraph_node` field in the streamed metadata:

```
# The "messages" stream mode streams LLM tokens with metadata
# Use version="v2" for a unified StreamPart format
for chunk in graph.stream(
 inputs,
 stream_mode="messages",
 version="v2",
):
 if chunk["type"] == "messages":
 msg, metadata = chunk["data"]
 # Filter the streamed tokens by the langgraph_node field in the metadata
 # to only include the tokens from the specified node
 if msg.content and metadata["langgraph_node"] == "some_node_name":
 ...

```

Extended example: streaming LLM tokens from specific nodes
```
from typing import TypedDict
from langgraph.graph import START, StateGraph
from langchain_openai import ChatOpenAI

model = ChatOpenAI(model="gpt-5.4-mini")

class State(TypedDict):
 topic: str
 joke: str
 poem: str

def write_joke(state: State):
 topic = state["topic"]
 joke_response = model.invoke(
 [{"role": "user", "content": f"Write a joke about {topic}"}]
 )
 return {"joke": joke_response.content}

def write_poem(state: State):
 topic = state["topic"]
 poem_response = model.invoke(
 [{"role": "user", "content": f"Write a short poem about {topic}"}]
 )
 return {"poem": poem_response.content}

graph = (
 StateGraph(State)
 .add_node(write_joke)
 .add_node(write_poem)
 # write both the joke and the poem concurrently
 .add_edge(START, "write_joke")
 .add_edge(START, "write_poem")
 .compile()
)

# The "messages" stream mode streams LLM tokens with metadata
# Use version="v2" for a unified StreamPart format
for chunk in graph.stream(
 {"topic": "cats"},
 stream_mode="messages",
 version="v2",
):
 if chunk["type"] == "messages":
 msg, metadata = chunk["data"]
 # Filter the streamed tokens by the langgraph_node field in the metadata
 # to only include the tokens from the write_poem node
 if msg.content and metadata["langgraph_node"] == "write_poem":
 print(msg.content, end="|", flush=True)

```

### [​](#custom-data)Custom data

To send **custom user-defined data** from inside a LangGraph node or tool, follow these steps:

{' ' * (self.list_depth - 1)}- Use [`get_stream_writer`](https://reference.langchain.com/python/langgraph/config/get_stream_writer) to access the stream writer and emit custom data.

{' ' * (self.list_depth - 1)}- Set `stream_mode="custom"` when calling `.stream()` or `.astream()` to get the custom data in the stream. You can combine multiple modes (e.g., `["updates", "custom"]`), but at least one must be `"custom"`.

**No [`get_stream_writer`](https://reference.langchain.com/python/langgraph/config/get_stream_writer) in async for Python < 3.11**
In async code running on Python < 3.11, [`get_stream_writer`](https://reference.langchain.com/python/langgraph/config/get_stream_writer) will not work.
Instead, add a `writer` parameter to your node or tool and pass it manually.
See [Async with Python < 3.11](#async) for usage examples.

{' ' * (self.list_depth - 1)}- node
{' ' * (self.list_depth - 1)}- tool
```
from typing import TypedDict
from langgraph.config import get_stream_writer
from langgraph.graph import StateGraph, START

class State(TypedDict):
 query: str
 answer: str

def node(state: State):
 # Get the stream writer to send custom data
 writer = get_stream_writer()
 # Emit a custom key-value pair (e.g., progress update)
 writer({"custom_key": "Generating custom data inside node"})
 return {"answer": "some data"}

graph = (
 StateGraph(State)
 .add_node(node)
 .add_edge(START, "node")
 .compile()
)

inputs = {"query": "example"}

# Set stream_mode="custom" to receive the custom data in the stream
for chunk in graph.stream(inputs, stream_mode="custom", version="v2"):
 if chunk["type"] == "custom":
 print(f"Custom event: {chunk['data']['custom_key']}")

```

```
from langchain.tools import tool
from langgraph.config import get_stream_writer

@tool
def query_database(query: str) -> str:
 """Query the database."""
 # Access the stream writer to send custom data
 writer = get_stream_writer()
 # Emit a custom key-value pair (e.g., progress update)
 writer({"data": "Retrieved 0/100 records", "type": "progress"})
 # perform query
 # Emit another custom key-value pair
 writer({"data": "Retrieved 100/100 records", "type": "progress"})
 return "some-answer"

graph = ... # define a graph that uses this tool

# Set stream_mode="custom" to receive the custom data in the stream
for chunk in graph.stream(inputs, stream_mode="custom", version="v2"):
 if chunk["type"] == "custom":
 print(f"{chunk['data']['type']}: {chunk['data']['data']}")

```

### [​](#subgraph-outputs)Subgraph outputs

To include outputs from [subgraphs](/oss/python/langgraph/use-subgraphs) in the streamed outputs, you can set `subgraphs=True` in the `.stream()` method of the parent graph. This will stream outputs from both the parent graph and any subgraphs.
The outputs will be streamed as tuples `(namespace, data)`, where `namespace` is a tuple with the path to the node where a subgraph is invoked, e.g. `("parent_node:<task_id>", "child_node:<task_id>")`.

{' ' * (self.list_depth - 1)}- v2 (LangGraph >= 1.1)
{' ' * (self.list_depth - 1)}- v1 (default)With `version="v2"`, subgraph events use the same `StreamPart` format. The `ns` field identifies the source:
```
for chunk in graph.stream(
 {"foo": "foo"},
 subgraphs=True,
 stream_mode="updates",
 version="v2",
):
 print(chunk["type"]) # "updates"
 print(chunk["ns"]) # () for root, ("node_name:<task_id>",) for subgraph
 print(chunk["data"]) # {"node_name": {"key": "value"}}

```

```
for chunk in graph.stream(
 {"foo": "foo"},
 # Set subgraphs=True to stream outputs from subgraphs
 subgraphs=True,
 stream_mode="updates",
):
 print(chunk)

```

Extended example: streaming from subgraphs
```
from langgraph.graph import START, StateGraph
from typing import TypedDict

# Define subgraph
class SubgraphState(TypedDict):
 foo: str # note that this key is shared with the parent graph state
 bar: str

def subgraph_node_1(state: SubgraphState):
 return {"bar": "bar"}

def subgraph_node_2(state: SubgraphState):
 return {"foo": state["foo"] + state["bar"]}

subgraph_builder = StateGraph(SubgraphState)
subgraph_builder.add_node(subgraph_node_1)
subgraph_builder.add_node(subgraph_node_2)
subgraph_builder.add_edge(START, "subgraph_node_1")
subgraph_builder.add_edge("subgraph_node_1", "subgraph_node_2")
subgraph = subgraph_builder.compile()

# Define parent graph
class ParentState(TypedDict):
 foo: str

def node_1(state: ParentState):
 return {"foo": "hi! " + state["foo"]}

builder = StateGraph(ParentState)
builder.add_node("node_1", node_1)
builder.add_node("node_2", subgraph)
builder.add_edge(START, "node_1")
builder.add_edge("node_1", "node_2")
graph = builder.compile()

for chunk in graph.stream(
 {"foo": "foo"},
 stream_mode="updates",
 # Set subgraphs=True to stream outputs from subgraphs
 subgraphs=True,
 version="v2",
):
 if chunk["type"] == "updates":
 if chunk["ns"]:
 print(f"Subgraph {chunk['ns']}: {chunk['data']}")
 else:
 print(f"Root: {chunk['data']}")

```

```
Root: {'node_1': {'foo': 'hi! foo'}}
Subgraph ('node_2:dfddc4ba-c3c5-6887-5012-a243b5b377c2',): {'subgraph_node_1': {'bar': 'bar'}}
Subgraph ('node_2:dfddc4ba-c3c5-6887-5012-a243b5b377c2',): {'subgraph_node_2': {'foo': 'hi! foobar'}}
Root: {'node_2': {'foo': 'hi! foobar'}}

```
**Note** that we are receiving not just the node updates, but we also the namespaces which tell us what graph (or subgraph) we are streaming from.

### [​](#checkpoints)Checkpoints

Use the `checkpoints` streaming mode to receive checkpoint events as the graph executes. Each checkpoint event has the same format as the output of `get_state()`. Requires a [checkpointer](/oss/python/langgraph/persistence).

```
from langgraph.checkpoint.memory import MemorySaver

graph = (
 StateGraph(State)
 .add_node(refine_topic)
 .add_node(generate_joke)
 .add_edge(START, "refine_topic")
 .add_edge("refine_topic", "generate_joke")
 .add_edge("generate_joke", END)
 .compile(checkpointer=MemorySaver())
)

config = {"configurable": {"thread_id": "1"}}

for chunk in graph.stream(
 {"topic": "ice cream"},
 config=config,
 stream_mode="checkpoints",
 version="v2",
):
 if chunk["type"] == "checkpoints":
 print(chunk["data"])

```

### [​](#tasks)Tasks

Use the `tasks` streaming mode to receive task start and finish events as the graph executes. Task events include information about which node is running, its results, and any errors. Requires a [checkpointer](/oss/python/langgraph/persistence).

```
from langgraph.checkpoint.memory import MemorySaver

graph = (
 StateGraph(State)
 .add_node(refine_topic)
 .add_node(generate_joke)
 .add_edge(START, "refine_topic")
 .add_edge("refine_topic", "generate_joke")
 .add_edge("generate_joke", END)
 .compile(checkpointer=MemorySaver())
)

config = {"configurable": {"thread_id": "1"}}

for chunk in graph.stream(
 {"topic": "ice cream"},
 config=config,
 stream_mode="tasks",
 version="v2",
):
 if chunk["type"] == "tasks":
 print(chunk["data"])

```

[]()

### [​](#debug)Debug

Use the `debug` streaming mode to stream as much information as possible throughout the execution of the graph. The streamed outputs include the name of the node as well as the full state.

```
for chunk in graph.stream(
 {"topic": "ice cream"},
 stream_mode="debug",
 version="v2",
):
 if chunk["type"] == "debug":
 print(chunk["data"])

```

The `debug` mode combines `checkpoints` and `tasks` events with additional metadata. Use `checkpoints` or `tasks` directly if you only need a subset of the debug information.

### [​](#multiple-modes-at-once)Multiple modes at once

You can pass a list as the `stream_mode` parameter to stream multiple modes at once.
With `version="v2"`, every chunk is a `StreamPart` dict. Use `chunk["type"]` to distinguish between modes:
v2v1
```
for chunk in graph.stream(inputs, stream_mode=["updates", "custom"], version="v2"):
 if chunk["type"] == "updates":
 for node_name, state in chunk["data"].items():
 print(f"Node `{node_name}` updated: {state}")
 elif chunk["type"] == "custom":
 print(f"Custom event: {chunk['data']}")

```

## [​](#advanced)Advanced

### [​](#use-with-any-llm)Use with any LLM

You can use `stream_mode="custom"` to stream data from **any LLM API**—even if that API does **not** implement the LangChain chat model interface.
This lets you integrate raw LLM clients or external services that provide their own streaming interfaces, making LangGraph highly flexible for custom setups.

```
from langgraph.config import get_stream_writer

def call_arbitrary_model(state):
 """Example node that calls an arbitrary model and streams the output"""
 # Get the stream writer to send custom data
 writer = get_stream_writer()
 # Assume you have a streaming client that yields chunks
 # Generate LLM tokens using your custom streaming client
 for chunk in your_custom_streaming_client(state["topic"]):
 # Use the writer to send custom data to the stream
 writer({"custom_llm_chunk": chunk})
 return {"result": "completed"}

graph = (
 StateGraph(State)
 .add_node(call_arbitrary_model)
 # Add other nodes and edges as needed
 .compile()
)
# Set stream_mode="custom" to receive the custom data in the stream
for chunk in graph.stream(
 {"topic": "cats"},
 stream_mode="custom",
 version="v2",
):
 if chunk["type"] == "custom":
 # The chunk data will contain the custom data streamed from the llm
 print(chunk["data"])

```

Extended example: streaming arbitrary chat model
```
import operator
import json

from typing import TypedDict
from typing_extensions import Annotated
from langgraph.graph import StateGraph, START

from openai import AsyncOpenAI

openai_client = AsyncOpenAI()
model_name = "gpt-5.4-mini"

async def stream_tokens(model_name: str, messages: list[dict]):
 response = await openai_client.chat.completions.create(
 messages=messages, model=model_name, stream=True
 )
 role = None
 async for chunk in response:
 delta = chunk.choices[0].delta

 if delta.role is not None:
 role = delta.role

 if delta.content:
 yield {"role": role, "content": delta.content}

# this is our tool
async def get_items(place: str) -> str:
 """Use this tool to list items one might find in a place you're asked about."""
 writer = get_stream_writer()
 response = ""
 async for msg_chunk in stream_tokens(
 model_name,
 [
 {
 "role": "user",
 "content": (
 "Can you tell me what kind of items "
 f"i might find in the following place: '{place}'. "
 "List at least 3 such items separating them by a comma. "
 "And include a brief description of each item."
 ),
 }
 ],
 ):
 response += msg_chunk["content"]
 writer(msg_chunk)

 return response

class State(TypedDict):
 messages: Annotated[list[dict], operator.add]

# this is the tool-calling graph node
async def call_tool(state: State):
 ai_message = state["messages"][-1]
 tool_call = ai_message["tool_calls"][-1]

 function_name = tool_call["function"]["name"]
 if function_name != "get_items":
 raise ValueError(f"Tool {function_name} not supported")

 function_arguments = tool_call["function"]["arguments"]
 arguments = json.loads(function_arguments)

 function_response = await get_items(**arguments)
 tool_message = {
 "tool_call_id": tool_call["id"],
 "role": "tool",
 "name": function_name,
 "content": function_response,
 }
 return {"messages": [tool_message]}

graph = (
 StateGraph(State)
 .add_node(call_tool)
 .add_edge(START, "call_tool")
 .compile()
)

```
Let’s invoke the graph with an [`AIMessage`](https://reference.langchain.com/python/langchain-core/messages/ai/AIMessage) that includes a tool call:
```
inputs = {
 "messages": [
 {
 "content": None,
 "role": "assistant",
 "tool_calls": [
 {
 "id": "1",
 "function": {
 "arguments": '{"place":"bedroom"}',
 "name": "get_items",
 },
 "type": "function",
 }
 ],
 }
 ]
}

async for chunk in graph.astream(
 inputs,
 stream_mode="custom",
 version="v2",
):
 if chunk["type"] == "custom":
 print(chunk["data"]["content"], end="|", flush=True)

```

### [​](#disable-streaming-for-specific-chat-models)Disable streaming for specific chat models

If your application mixes models that support streaming with those that do not, you may need to explicitly disable streaming for
models that do not support it.
Set `streaming=False` when initializing the model.

{' ' * (self.list_depth - 1)}- init_chat_model
{' ' * (self.list_depth - 1)}- Chat model interface
```
from langchain.chat_models import init_chat_model

model = init_chat_model(
 "claude-sonnet-4-6",
 # Set streaming=False to disable streaming for the chat model
 streaming=False
)

```

```
from langchain_openai import ChatOpenAI

# Set streaming=False to disable streaming for the chat model
model = ChatOpenAI(model="o1-preview", streaming=False)

```

Not all chat model integrations support the `streaming` parameter. If your model doesn’t support it, use `disable_streaming=True` instead. This parameter is available on all chat models via the base class.

### [​](#migrate-to-v2)Migrate to v2

The v2 streaming format (used throughout this page) provides a unified output format. Here’s a summary of the key differences and how to migrate:

| 
 | Scenario | v1 (default) | v2 (`version="v2"`)
 | Single stream mode | Raw data (dict) | `StreamPart` dict with `type`, `ns`, `data`
 | Multiple stream modes | `(mode, data)` tuples | Same `StreamPart` dict, filter on `chunk["type"]`
 | Subgraph streaming | `(namespace, data)` tuples | Same `StreamPart` dict, check `chunk["ns"]`
 | Multiple modes + subgraphs | `(namespace, mode, data)` triples | Same `StreamPart` dict
 | `invoke()` return type | Plain dict (state) | `GraphOutput` with `.value` and `.interrupts`
 | Interrupt location (stream) | `__interrupt__` key in state dict | `interrupts` field on `values` stream parts
 | Interrupt location (invoke) | `__interrupt__` key in result dict | `.interrupts` attribute on `GraphOutput`
 | Pydantic/dataclass output | Returns plain dict | Coerces to model/dataclass instance

#### [​](#v2-invoke-format)v2 invoke format

When you pass `version="v2"` to `invoke()` or `ainvoke()`, it returns a [`GraphOutput`](https://reference.langchain.com/python/langgraph/types/GraphOutput) object with `.value` and `.interrupts` attributes:

```
from langgraph.types import GraphOutput

result = graph.invoke(inputs, version="v2")

assert isinstance(result, GraphOutput)
result.value # your output — dict, Pydantic model, or dataclass
result.interrupts # tuple[Interrupt, ...], empty if none occurred

```

With any stream mode other than the default `"values"`, `invoke(..., stream_mode="updates", version="v2")` returns `list[StreamPart]` instead of `list[tuple]`.
Dict-style access on `GraphOutput` (`result["key"]`, `"key" in result`, `result["__interrupt__"]`) still works for backwards compatibility but is **deprecated** and will be removed in a future version. Migrate to `result.value` and `result.interrupts`.
This separates state from interrupt metadata. With v1, interrupts are embedded in the returned dict under `__interrupt__`:
v2 (new)v1 (current default)
```
config = {"configurable": {"thread_id": "thread-1"}}
result = graph.invoke(inputs, config=config, version="v2")

if result.interrupts:
 print(result.interrupts[0].value)
 graph.invoke(Command(resume=True), config=config, version="v2")

```

#### [​](#pydantic-and-dataclass-state-coercion)Pydantic and dataclass state coercion

When your graph state is a Pydantic model or dataclass, v2 `values` mode automatically coerces output to the correct type:

```
from pydantic import BaseModel
from typing import Annotated
import operator

class MyState(BaseModel):
 value: str
 items: Annotated[list[str], operator.add]

# With version="v2", chunk["data"] is a MyState instance
for chunk in graph.stream(
 {"value": "x", "items": []}, stream_mode="values", version="v2"
):
 print(type(chunk["data"])) # <class 'MyState'>

```

[]()

### [​](#async-with-python-<-3-11)Async with Python < 3.11

In Python versions < 3.11, [asyncio tasks](https://docs.python.org/3/library/asyncio-task.html#asyncio.create_task) do not support the `context` parameter.
This limits LangGraph ability to automatically propagate context, and affects LangGraph’s streaming mechanisms in two key ways:

{' ' * (self.list_depth - 1)}- You **must** explicitly pass [`RunnableConfig`](https://python.langchain.com/docs/concepts/runnables/#runnableconfig) into async LLM calls (e.g., `ainvoke()`), as callbacks are not automatically propagated.

{' ' * (self.list_depth - 1)}- You **cannot** use [`get_stream_writer`](https://reference.langchain.com/python/langgraph/config/get_stream_writer) in async nodes or tools—you must pass a `writer` argument directly.

Extended example: async LLM call with manual config
```
from typing import TypedDict
from langgraph.graph import START, StateGraph
from langchain.chat_models import init_chat_model

model = init_chat_model(model="gpt-5.4-mini")

class State(TypedDict):
 topic: str
 joke: str

# Accept config as an argument in the async node function
async def call_model(state, config):
 topic = state["topic"]
 print("Generating joke...")
 # Pass config to model.ainvoke() to ensure proper context propagation
 joke_response = await model.ainvoke(
 [{"role": "user", "content": f"Write a joke about {topic}"}],
 config,
 )
 return {"joke": joke_response.content}

graph = (
 StateGraph(State)
 .add_node(call_model)
 .add_edge(START, "call_model")
 .compile()
)

# Set stream_mode="messages" to stream LLM tokens
async for chunk in graph.astream(
 {"topic": "ice cream"},
 stream_mode="messages",
 version="v2",
):
 if chunk["type"] == "messages":
 message_chunk, metadata = chunk["data"]
 if message_chunk.content:
 print(message_chunk.content, end="|", flush=True)

```

Extended example: async custom streaming with stream writer
```
from typing import TypedDict
from langgraph.types import StreamWriter

class State(TypedDict):
 topic: str
 joke: str

# Add writer as an argument in the function signature of the async node or tool
# LangGraph will automatically pass the stream writer to the function
async def generate_joke(state: State, writer: StreamWriter):
 writer({"custom_key": "Streaming custom data while generating a joke"})
 return {"joke": f"This is a joke about {state['topic']}"}

graph = (
 StateGraph(State)
 .add_node(generate_joke)
 .add_edge(START, "generate_joke")
 .compile()
)

# Set stream_mode="custom" to receive the custom data in the stream #
async for chunk in graph.astream(
 {"topic": "ice cream"},
 stream_mode="custom",
 version="v2",
):
 if chunk["type"] == "custom":
 print(chunk["data"])

```

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langgraph/streaming.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Event streamingPrevious](/oss/python/langgraph/event-streaming)[InterruptsNext](/oss/python/langgraph/interrupts)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
