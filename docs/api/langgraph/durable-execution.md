# Durable execution

[Capabilities](/oss/python/langgraph/persistence)

# Durable execution
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.**Durable execution** is a technique in which a process or workflow saves its progress at key points, allowing it to pause and later resume exactly where it left off. This is particularly useful in scenarios that require [human-in-the-loop](/oss/python/langgraph/interrupts), where users can inspect, validate, or modify the process before continuing, and in long-running tasks that might encounter interruptions or errors (e.g., calls to an LLM timing out). By preserving completed work, durable execution enables a process to resume without reprocessing previous steps — even after a significant delay (e.g., a week later).
LangGraph’s built-in [persistence](/oss/python/langgraph/persistence) layer provides durable execution for workflows, ensuring that the state of each execution step is saved to a durable store. This capability guarantees that if a workflow is interrupted — whether by a system failure or for [human-in-the-loop](/oss/python/langgraph/interrupts) interactions — it can be resumed from its last recorded state.
If you are using LangGraph with a checkpointer, you already have durable execution enabled. You can pause and resume workflows at any point, even after interruptions or failures.
To make the most of durable execution, ensure that your workflow is designed to be [deterministic](#determinism-and-consistent-replay) and [idempotent](#determinism-and-consistent-replay) and wrap any side effects or non-deterministic operations inside [tasks](/oss/python/langgraph/functional-api#task). You can use [tasks](/oss/python/langgraph/functional-api#task) from both the [StateGraph (Graph API)](/oss/python/langgraph/graph-api) and the [Functional API](/oss/python/langgraph/functional-api).

## [​](#requirements)Requirements

To leverage durable execution in LangGraph, you need to:

{' ' * (self.list_depth - 1)}- 
Enable [persistence](/oss/python/langgraph/persistence) in your workflow by specifying a [checkpointer](/oss/python/langgraph/persistence#checkpointer-libraries) that will save workflow progress.

{' ' * (self.list_depth - 1)}- 
Specify a [thread identifier](/oss/python/langgraph/persistence#threads) when executing a workflow. This will track the execution history for a particular instance of the workflow.

{' ' * (self.list_depth - 1)}- 
Wrap any non-deterministic operations (e.g., random number generation) or operations with side effects (e.g., file writes, API calls) inside [tasks](https://reference.langchain.com/python/langgraph/func/task) to ensure that when a workflow is resumed, these operations are not repeated for the particular run, and instead their results are retrieved from the persistence layer. For more information, see [Determinism and Consistent Replay](#determinism-and-consistent-replay).

## [​](#determinism-and-consistent-replay)Determinism and consistent replay

When you resume a workflow run, the code does **NOT** resume from the **same line of code** where execution stopped; instead, it will identify an appropriate [starting point](#starting-points-for-resuming-workflows) from which to pick up where it left off. This means that the workflow will replay all steps from the [starting point](#starting-points-for-resuming-workflows) until it reaches the point where it was stopped.
As a result, when you are writing a workflow for durable execution, you must wrap any non-deterministic operations (e.g., random number generation) and any operations with side effects (e.g., file writes, API calls) inside [tasks](/oss/python/langgraph/functional-api#task) or [nodes](/oss/python/langgraph/graph-api#nodes).
To ensure that your workflow is deterministic and can be consistently replayed, follow these guidelines:

{' ' * (self.list_depth - 1)}- **Avoid Repeating Work**: If a [node](/oss/python/langgraph/graph-api#nodes) contains multiple operations with side effects (e.g., logging, file writes, or network calls), wrap each operation in a separate **task**. This ensures that when the workflow is resumed, the operations are not repeated, and their results are retrieved from the persistence layer.

{' ' * (self.list_depth - 1)}- **Encapsulate Non-Deterministic Operations:** Wrap any code that might yield non-deterministic results (e.g., random number generation) inside **tasks** or **nodes**. This ensures that, upon resumption, the workflow follows the exact recorded sequence of steps with the same outcomes.

{' ' * (self.list_depth - 1)}- **Use Idempotent Operations**: When possible ensure that side effects (e.g., API calls, file writes) are idempotent. This means that if an operation is retried after a failure in the workflow, it will have the same effect as the first time it was executed. This is particularly important for operations that result in data writes. In the event that a **task** starts but fails to complete successfully, the workflow’s resumption will re-run the **task**, relying on recorded outcomes to maintain consistency. Use idempotency keys or verify existing results to avoid unintended duplication, ensuring a smooth and predictable workflow execution.

For some examples of pitfalls to avoid, see the [Common Pitfalls](/oss/python/langgraph/functional-api#common-pitfalls) section in the functional API, which shows
how to structure your code using **tasks** to avoid these issues. The same principles apply to the [StateGraph (Graph API)](https://reference.langchain.com/python/langgraph/graph/state/StateGraph).

## [​](#durability-modes)Durability modes

LangGraph supports three durability modes that allow you to balance performance and data consistency based on your application’s requirements. A higher durability mode adds more overhead to the workflow execution. You can specify the durability mode when calling any graph execution method:

```
graph.stream(
 {"input": "test"},
 durability="sync"
)

```

The durability modes, from least to most durable, are as follows:

{' ' * (self.list_depth - 1)}- `"exit"`: LangGraph persists changes only when graph execution exits either successfully, with an error, or due to a human in the loop interrupt. This provides the best performance for long-running graphs but means intermediate state is not saved, so you cannot recover from system failures (like process crashes) that occur mid-execution.

{' ' * (self.list_depth - 1)}- `"async"`: LangGraph persists changes asynchronously while the next step executes. This provides good performance and durability, but there’s a small risk that LangGraph does not write checkpoints if the process crashes during execution.

{' ' * (self.list_depth - 1)}- `"sync"`: LangGraph persists changes synchronously before the next step starts. This ensures that LangGraph writes every checkpoint before continuing execution, providing high durability at the cost of some performance overhead.

## [​](#using-tasks-in-nodes)Using tasks in nodes

If a [node](/oss/python/langgraph/graph-api#nodes) contains multiple operations, you may find it easier to convert each operation into a **task** rather than refactor the operations into individual nodes.

{' ' * (self.list_depth - 1)}- Original
{' ' * (self.list_depth - 1)}- With task
```
from typing import NotRequired
from typing_extensions import TypedDict
from langchain_core.utils.uuid import uuid7

from langgraph.checkpoint.memory import InMemorySaver
from langgraph.graph import StateGraph, START, END
import requests

# Define a TypedDict to represent the state
class State(TypedDict):
 url: str
 result: NotRequired[str]

def call_api(state: State):
 """Example node that makes an API request."""
 result = requests.get(state['url']).text[:100] # Side-effect #
 return {
 "result": result
 }

# Create a StateGraph builder and add a node for the call_api function
builder = StateGraph(State)
builder.add_node("call_api", call_api)

# Connect the start and end nodes to the call_api node
builder.add_edge(START, "call_api")
builder.add_edge("call_api", END)

# Specify a checkpointer
checkpointer = InMemorySaver()

# Compile the graph with the checkpointer
graph = builder.compile(checkpointer=checkpointer)

# Define a config with a thread ID.
thread_id = str(uuid7())
config = {"configurable": {"thread_id": thread_id}}

# Invoke the graph
graph.invoke({"url": "https://www.example.com"}, config)

```

```
from typing import NotRequired
from typing_extensions import TypedDict
from langchain_core.utils.uuid import uuid7

from langgraph.checkpoint.memory import InMemorySaver
from langgraph.func import task
from langgraph.graph import StateGraph, START, END
import requests

# Define a TypedDict to represent the state
class State(TypedDict):
 urls: list[str]
 result: NotRequired[list[str]]

@task
def _make_request(url: str):
 """Make a request."""
 return requests.get(url).text[:100]

def call_api(state: State):
 """Example node that makes an API request."""
 requests = [_make_request(url) for url in state['urls']]
 results = [request.result() for request in requests]
 return {
 "results": results
 }

# Create a StateGraph builder and add a node for the call_api function
builder = StateGraph(State)
builder.add_node("call_api", call_api)

# Connect the start and end nodes to the call_api node
builder.add_edge(START, "call_api")
builder.add_edge("call_api", END)

# Specify a checkpointer
checkpointer = InMemorySaver()

# Compile the graph with the checkpointer
graph = builder.compile(checkpointer=checkpointer)

# Define a config with a thread ID.
thread_id = str(uuid7())
config = {"configurable": {"thread_id": thread_id}}

# Invoke the graph
graph.invoke({"urls": ["https://www.example.com"]}, config)

```

## [​](#resuming-workflows)Resuming workflows

Once you have enabled durable execution in your workflow, you can resume execution for the following scenarios:

{' ' * (self.list_depth - 1)}- **Pausing and Resuming Workflows:** Use the [interrupt](https://reference.langchain.com/python/langgraph/types/interrupt) function to pause a workflow at specific points and the [`Command`](https://reference.langchain.com/python/langgraph/types/Command) primitive to resume it with updated state. See [**Interrupts**](/oss/python/langgraph/interrupts) for more details.

{' ' * (self.list_depth - 1)}- **Recovering from Failures:** Automatically resume workflows from the last successful checkpoint after an exception (e.g., LLM provider outage). This involves executing the workflow with the same thread identifier by providing it with a `None` as the input value (see this [example](/oss/python/langgraph/use-functional-api#resuming-after-an-error) with the functional API).

## [​](#starting-points-for-resuming-workflows)Starting points for resuming workflows

{' ' * (self.list_depth - 1)}- If you’re using a [StateGraph (Graph API)](/oss/python/langgraph/graph-api), the starting point is the beginning of the [**node**](/oss/python/langgraph/graph-api#nodes) where execution stopped.

{' ' * (self.list_depth - 1)}- If you’re making a subgraph call inside a node, the starting point will be the **parent** node that called the subgraph that was halted.
Inside the subgraph, the starting point will be the specific [**node**](/oss/python/langgraph/graph-api#nodes) where execution stopped.

{' ' * (self.list_depth - 1)}- If you’re using the Functional API, the starting point is the beginning of the [**entrypoint**](/oss/python/langgraph/functional-api#entrypoint) where execution stopped.

## [​](#graceful-shutdown)Graceful shutdown

Requires `langgraph>=1.2`.
Graceful shutdown lets you stop an in-flight graph run cooperatively—after the current superstep completes—and save a resumable checkpoint. This is useful for handling SIGTERM signals or any external supervisor that needs to reclaim resources without losing work.
Create a [`RunControl`](https://reference.langchain.com/python/langgraph/runtime/RunControl) and pass it as `control=` to `invoke` or `stream`. Call `request_drain()` from any thread to signal that the run should stop:

```
from langgraph.runtime import RunControl
from langgraph.errors import GraphDrained

control = RunControl()

# In a signal handler or supervisor:
# control.request_drain("sigterm")

try:
 result = graph.invoke(inputs, config, control=control)
except GraphDrained as e:
 # The graph stopped early and saved a checkpoint.
 # Resume later with the same config.
 print(f"Drained: {e.reason}")

```

### [​](#semantics)Semantics

Drain is cooperative and operates between supersteps, never preempting work that is already running:

| 
 | Scenario | Behavior
 | Node mid-execution | Runs to completion. Drain takes effect on the next superstep.
 | Node with a retry policy currently retrying | Retry loop runs to exhaustion or success. Drain takes effect after.
 | Graph finishes naturally on the same tick as drain | Returns normally. Inspect `control.drain_requested` to distinguish from a normal run.
 | More supersteps remain | Raises `GraphDrained(reason)`. Checkpoint is saved and resumable.
 | Subgraph requests drain | `GraphDrained` bubbles up through the parent and stops it at its own next superstep boundary.

### [​](#resume-after-drain)Resume after drain

Resume a drained run with `invoke(None, config)` using the same `thread_id`:

```
result = graph.invoke(None, config)

```

### [​](#read-drain-state-inside-a-node)Read drain state inside a node

Access drain state through the `runtime` parameter to adjust node behavior before the superstep boundary is reached:

```
from langgraph.runtime import Runtime

async def my_node(state: State, runtime: Runtime) -> State:
 if runtime.drain_requested:
 # Skip expensive work and return a minimal result
 return {"status": "skipped", "reason": runtime.drain_reason}
 return {"status": await do_work()}

```

### [​](#sigterm-hook-pattern)SIGTERM hook pattern

The recommended pattern for handling process shutdown:

```
import signal
from langgraph.runtime import RunControl
from langgraph.errors import GraphDrained

control = RunControl()
signal.signal(signal.SIGTERM, lambda *_: control.request_drain("sigterm"))

try:
 result = graph.invoke(inputs, config, control=control)
except GraphDrained as e:
 log.info("graph drained: %s", e.reason)
 # Resume on next startup with the same config

```

`request_drain()` does not cancel running asyncio tasks or kill threads. For a hard upper bound, pair drain with a graceful timeout and task cancellation.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langgraph/durable-execution.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[PersistencePrevious](/oss/python/langgraph/persistence)[Fault toleranceNext](/oss/python/langgraph/fault-tolerance)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
