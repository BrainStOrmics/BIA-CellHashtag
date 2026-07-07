# Functional API overview

[LangGraph APIs](/oss/python/langgraph/choosing-apis)[Functional API](/oss/python/langgraph/functional-api)

# Functional API overview
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.The **Functional API** allows you to add LangGraph’s key features ([persistence](/oss/python/langgraph/persistence), [memory](/oss/python/langgraph/add-memory), [human-in-the-loop](/oss/python/langgraph/interrupts), and [streaming](/oss/python/langgraph/streaming)) to your applications with minimal changes to your existing code.
It is designed to integrate these features into existing code that may use standard language primitives for branching and control flow, such as `if` statements, `for` loops, and function calls. Unlike many data orchestration frameworks that require restructuring code into an explicit pipeline or DAG, the Functional API allows you to incorporate these capabilities without enforcing a rigid execution model.
The Functional API uses two key building blocks:

{' ' * (self.list_depth - 1)}- **`@entrypoint`**: Marks a function as the starting point of a workflow, encapsulating logic and managing execution flow, including handling long-running tasks and interrupts.

{' ' * (self.list_depth - 1)}- **[`@task`](https://reference.langchain.com/python/langgraph/func/task)**: Represents a discrete unit of work, such as an API call or data processing step, that can be executed asynchronously within an entrypoint. Tasks return a future-like object that can be awaited or resolved synchronously.

This provides a minimal abstraction for building workflows with state management and streaming.
For information on how to use the functional API, see [Use Functional API](/oss/python/langgraph/use-functional-api).

## [​](#functional-api-vs-graph-api)Functional API vs. Graph API

For users who prefer a more declarative approach, LangGraph’s [Graph API](/oss/python/langgraph/graph-api) allows you to define workflows using a Graph paradigm. Both APIs share the same underlying runtime, so you can use them together in the same application.
Here are some key differences:

{' ' * (self.list_depth - 1)}- **Control flow**: The Functional API does not require thinking about graph structure. You can use standard Python constructs to define workflows. This will usually trim the amount of code you need to write.

{' ' * (self.list_depth - 1)}- **Short-term memory**: The **GraphAPI** requires declaring a [**State**](/oss/python/langgraph/graph-api#state) and may require defining [**reducers**](/oss/python/langgraph/graph-api#reducers) to manage updates to the graph state. `@entrypoint` and `@tasks` do not require explicit state management as their state is scoped to the function and is not shared across functions.

{' ' * (self.list_depth - 1)}- **Checkpointing**: Both APIs generate and use checkpoints. In the **Graph API** a new checkpoint is generated after every [superstep](/oss/python/langgraph/graph-api). In the **Functional API**, when tasks are executed, their results are saved to an existing checkpoint associated with the given entrypoint instead of creating a new checkpoint.

{' ' * (self.list_depth - 1)}- **Visualization**: The Graph API makes it easy to visualize the workflow as a graph which can be useful for debugging, understanding the workflow, and sharing with others. The Functional API does not support visualization as the graph is dynamically generated during runtime.

## [​](#example)Example

Below we demonstrate a simple application that writes an essay and [interrupts](/oss/python/langgraph/interrupts) to request human review.

```
from langgraph.checkpoint.memory import InMemorySaver
from langgraph.func import entrypoint, task
from langgraph.types import interrupt

@task
def write_essay(topic: str) -> str:
 """Write an essay about the given topic."""
 time.sleep(1) # A placeholder for a long-running task.
 return f"An essay about topic: {topic}"

@entrypoint(checkpointer=InMemorySaver())
def workflow(topic: str) -> dict:
 """A simple workflow that writes an essay and asks for a review."""
 essay = write_essay("cat").result()
 is_approved = interrupt({
 # Any json-serializable payload provided to interrupt as argument.
 # It will be surfaced on the client side as an Interrupt when streaming data
 # from the workflow.
 "essay": essay, # The essay we want reviewed.
 # We can add any additional information that we need.
 # For example, introduce a key called "action" with some instructions.
 "action": "Please approve/reject the essay",
 })

 return {
 "essay": essay, # The essay that was generated
 "is_approved": is_approved, # Response from HIL
 }

```

Detailed ExplanationThis workflow will write an essay about the topic “cat” and then pause to get a review from a human. The workflow can be interrupted for an indefinite amount of time until a review is provided.When the workflow is resumed, it executes from the very start, but because the result of the `writeEssay` task was already saved, the task result will be loaded from the checkpoint instead of being recomputed.
```
import time
from langchain_core.utils.uuid import uuid7
from langgraph.func import entrypoint, task
from langgraph.types import interrupt
from langgraph.checkpoint.memory import InMemorySaver

@task
def write_essay(topic: str) -> str:
 """Write an essay about the given topic."""
 time.sleep(1) # This is a placeholder for a long-running task.
 return f"An essay about topic: {topic}"

@entrypoint(checkpointer=InMemorySaver())
def workflow(topic: str) -> dict:
 """A simple workflow that writes an essay and asks for a review."""
 essay = write_essay("cat").result()
 is_approved = interrupt(
 {
 # Any json-serializable payload provided to interrupt as argument.
 # It will be surfaced on the client side as an Interrupt when streaming data
 # from the workflow.
 "essay": essay, # The essay we want reviewed.
 # We can add any additional information that we need.
 # For example, introduce a key called "action" with some instructions.
 "action": "Please approve/reject the essay",
 }
 )
 return {
 "essay": essay, # The essay that was generated
 "is_approved": is_approved, # Response from HIL
 }

thread_id = str(uuid7())
config = {"configurable": {"thread_id": thread_id}}
for item in workflow.stream("cat", config):
 print(item)
# > {'write_essay': 'An essay about topic: cat'}
# > {
# > '__interrupt__': (
# > Interrupt(
# > value={
# > 'essay': 'An essay about topic: cat',
# > 'action': 'Please approve/reject the essay'
# > },
# > id='b9b2b9d788f482663ced6dc755c9e981'
# > ),
# > )
# > }

```
An essay has been written and is ready for review. Once the review is provided, we can resume the workflow:
```
from langgraph.types import Command

# Get review from a user (e.g., via a UI)
# In this case, we're using a bool, but this can be any json-serializable value.
human_review = True

for item in workflow.stream(Command(resume=human_review), config):
 print(item)

```

```
{'workflow': {'essay': 'An essay about topic: cat', 'is_approved': False}}

```
The workflow has been completed and the review has been added to the essay.

## [​](#entrypoint)Entrypoint

The [`@entrypoint`](https://reference.langchain.com/python/langgraph/func/entrypoint) decorator can be used to create a workflow from a function. It encapsulates workflow logic and manages execution flow, including handling *long-running tasks* and [interrupts](/oss/python/langgraph/interrupts).

### [​](#definition)Definition

An **entrypoint** is defined by decorating a function with the `@entrypoint` decorator.
The function **must accept a single positional argument**, which serves as the workflow input. If you need to pass multiple pieces of data, use a dictionary as the input type for the first argument.
Decorating a function with an `entrypoint` produces a [`Pregel`](https://reference.langchain.com/python/langgraph/pregel/#langgraph.pregel.Pregel.stream) instance which helps to manage the execution of the workflow (e.g., handles streaming, resumption, and checkpointing).
You will usually want to pass a **checkpointer** to the `@entrypoint` decorator to enable persistence and use features like **human-in-the-loop**.

{' ' * (self.list_depth - 1)}- Sync
{' ' * (self.list_depth - 1)}- Async
```
from langgraph.func import entrypoint

@entrypoint(checkpointer=checkpointer)
def my_workflow(some_input: dict) -> int:
 # some logic that may involve long-running tasks like API calls,
 # and may be interrupted for human-in-the-loop.
 ...
 return result

```

```
from langgraph.func import entrypoint

@entrypoint(checkpointer=checkpointer)
async def my_workflow(some_input: dict) -> int:
 # some logic that may involve long-running tasks like API calls,
 # and may be interrupted for human-in-the-loop
 ...
 return result

```

**Serialization**
The **inputs** and **outputs** of entrypoints must be JSON-serializable to support checkpointing. Please see the [serialization](#serialization) section for more details.

### [​](#injectable-parameters)Injectable parameters

When declaring an `entrypoint`, you can request access to additional parameters that will be injected automatically at runtime. These parameters include:

| 
 | Parameter | Description
 | **previous** | Access the state associated with the previous `checkpoint` for the given thread. See [short-term-memory](#short-term-memory).
 | **store** | An instance of [BaseStore][langgraph.store.base.BaseStore]. Useful for [long-term memory](/oss/python/langgraph/use-functional-api#long-term-memory).
 | **writer** | Use to access the StreamWriter when working with Async Python < 3.11. See [streaming with functional API for details](/oss/python/langgraph/use-functional-api#streaming).
 | **config** | For accessing run time configuration. See [RunnableConfig](https://python.langchain.com/docs/concepts/runnables/#runnableconfig) for information.
Declare the parameters with the appropriate name and type annotation.

Requesting Injectable Parameters
```
from langchain_core.runnables import RunnableConfig
from langgraph.func import entrypoint
from langgraph.store.base import BaseStore
from langgraph.store.memory import InMemoryStore
from langgraph.checkpoint.memory import InMemorySaver
from langgraph.types import StreamWriter

in_memory_checkpointer = InMemorySaver(...)
in_memory_store = InMemoryStore(...) # An instance of InMemoryStore for long-term memory

@entrypoint(
 checkpointer=in_memory_checkpointer, # Specify the checkpointer
 store=in_memory_store # Specify the store
)
def my_workflow(
 some_input: dict, # The input (e.g., passed via `invoke`)
 *,
 previous: Any = None, # For short-term memory
 store: BaseStore, # For long-term memory
 writer: StreamWriter, # For streaming custom data
 config: RunnableConfig # For accessing the configuration passed to the entrypoint
) -> ...:

```

### [​](#executing)Executing

Using the [`@entrypoint`](#entrypoint) yields a [`Pregel`](https://reference.langchain.com/python/langgraph/pregel/#langgraph.pregel.Pregel.stream) object that can be executed using the `invoke`, `ainvoke`, `stream`, and `astream` methods.

{' ' * (self.list_depth - 1)}- Invoke
{' ' * (self.list_depth - 1)}- Async Invoke
{' ' * (self.list_depth - 1)}- Stream
{' ' * (self.list_depth - 1)}- Async Stream
```
config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}
my_workflow.invoke(some_input, config) # Wait for the result synchronously

```

```
config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}
await my_workflow.ainvoke(some_input, config) # Await result asynchronously

```

```
config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}

for chunk in my_workflow.stream(some_input, config):
 print(chunk)

```

```
config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}

async for chunk in my_workflow.astream(some_input, config):
 print(chunk)

```

### [​](#resuming)Resuming

Resuming an execution after an [interrupt](https://reference.langchain.com/python/langgraph/types/interrupt) can be done by passing a **resume** value to the [`Command`](https://reference.langchain.com/python/langgraph/types/Command) primitive.

{' ' * (self.list_depth - 1)}- Invoke
{' ' * (self.list_depth - 1)}- Async Invoke
{' ' * (self.list_depth - 1)}- Stream
{' ' * (self.list_depth - 1)}- Async Stream
```
from langgraph.types import Command

config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}

my_workflow.invoke(Command(resume=some_resume_value), config)

```

```
from langgraph.types import Command

config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}

await my_workflow.ainvoke(Command(resume=some_resume_value), config)

```

```
from langgraph.types import Command

config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}

for chunk in my_workflow.stream(Command(resume=some_resume_value), config):
 print(chunk)

```

```
from langgraph.types import Command

config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}

async for chunk in my_workflow.astream(Command(resume=some_resume_value), config):
 print(chunk)

```

**Resuming after an error**
To resume after an error, run the `entrypoint` with a `None` and the same **thread id** (config).
This assumes that the underlying **error** has been resolved and execution can proceed successfully.

{' ' * (self.list_depth - 1)}- Invoke
{' ' * (self.list_depth - 1)}- Async Invoke
{' ' * (self.list_depth - 1)}- Stream
{' ' * (self.list_depth - 1)}- Async Stream
```

config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}

my_workflow.invoke(None, config)

```

```

config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}

await my_workflow.ainvoke(None, config)

```

```

config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}

for chunk in my_workflow.stream(None, config):
 print(chunk)

```

```

config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}

async for chunk in my_workflow.astream(None, config):
 print(chunk)

```

### [​](#short-term-memory)Short-term memory

When an `entrypoint` is defined with a `checkpointer`, it stores information between successive invocations on the same **thread id** in [checkpoints](/oss/python/langgraph/persistence#checkpoints).
This allows accessing the state from the previous invocation using the `previous` parameter.
By default, the `previous` parameter is the return value of the previous invocation.

```
@entrypoint(checkpointer=checkpointer)
def my_workflow(number: int, *, previous: Any = None) -> int:
 previous = previous or 0
 return number + previous

config = {
 "configurable": {
 "thread_id": "some_thread_id"
 }
}

my_workflow.invoke(1, config) # 1 (previous was None)
my_workflow.invoke(2, config) # 3 (previous was 1 from the previous invocation)

```

#### [​](#entrypoint-final)`entrypoint.final`

[`entrypoint.final`](https://reference.langchain.com/python/langgraph/func/entrypoint/final) is a special primitive that can be returned from an entrypoint and allows **decoupling** the value that is **saved in the checkpoint** from the **return value of the entrypoint**.
The first value is the return value of the entrypoint, and the second value is the value that will be saved in the checkpoint. The type annotation is `entrypoint.final[return_type, save_type]`.

```
@entrypoint(checkpointer=checkpointer)
def my_workflow(number: int, *, previous: Any = None) -> entrypoint.final[int, int]:
 previous = previous or 0
 # This will return the previous value to the caller, saving
 # 2 * number to the checkpoint, which will be used in the next invocation
 # for the `previous` parameter.
 return entrypoint.final(value=previous, save=2 * number)

config = {
 "configurable": {
 "thread_id": "1"
 }
}

my_workflow.invoke(3, config) # 0 (previous was None)
my_workflow.invoke(1, config) # 6 (previous was 3 * 2 from the previous invocation)

```

## [​](#task)Task

A **task** represents a discrete unit of work, such as an API call or data processing step. It has two key characteristics:

{' ' * (self.list_depth - 1)}- **Asynchronous Execution**: Tasks are designed to be executed asynchronously, allowing multiple operations to run concurrently without blocking.

{' ' * (self.list_depth - 1)}- **Checkpointing**: Task results are saved to a checkpoint, enabling resumption of the workflow from the last saved state. (See [persistence](/oss/python/langgraph/persistence) for more details).

### [​](#definition-2)Definition

Tasks are defined using the `@task` decorator, which wraps a regular Python function.

```
from langgraph.func import task

@task()
def slow_computation(input_value):
 # Simulate a long-running operation
 ...
 return result

```

**Serialization**
The **outputs** of tasks must be JSON-serializable to support checkpointing.

### [​](#execution)Execution

**Tasks** can only be called from within an **entrypoint**, another **task**, or a [state graph node](/oss/python/langgraph/graph-api#nodes).
Tasks *cannot* be called directly from the main application code.
When you call a **task**, it returns *immediately* with a future object. A future is a placeholder for a result that will be available later.
To obtain the result of a **task**, you can either wait for it synchronously (using `result()`) or await it asynchronously (using `await`).

{' ' * (self.list_depth - 1)}- Synchronous Invocation
{' ' * (self.list_depth - 1)}- Asynchronous Invocation
```
@entrypoint(checkpointer=checkpointer)
def my_workflow(some_input: int) -> int:
 future = slow_computation(some_input)
 return future.result() # Wait for the result synchronously

```

```
@entrypoint(checkpointer=checkpointer)
async def my_workflow(some_input: int) -> int:
 return await slow_computation(some_input) # Await result asynchronously

```

## [​](#when-to-use-a-task)When to use a task

**Tasks** are useful in the following scenarios:

{' ' * (self.list_depth - 1)}- **Checkpointing**: When you need to save the result of a long-running operation to a checkpoint, so you don’t need to recompute it when resuming the workflow.

{' ' * (self.list_depth - 1)}- **Human-in-the-loop**: If you’re building a workflow that requires human intervention, you MUST use **tasks** to encapsulate any randomness (e.g., API calls) to ensure that the workflow can be resumed correctly. See the [determinism](#determinism) section for more details.

{' ' * (self.list_depth - 1)}- **Parallel Execution**: For I/O-bound tasks, **tasks** enable parallel execution, allowing multiple operations to run concurrently without blocking (e.g., calling multiple APIs).

{' ' * (self.list_depth - 1)}- **Observability**: Wrapping operations in **tasks** provides a way to track the progress of the workflow and monitor the execution of individual operations using [LangSmith](/langsmith/home).

{' ' * (self.list_depth - 1)}- **Retryable Work**: When work needs to be retried to handle failures or inconsistencies, **tasks** provide a way to encapsulate and manage the retry logic.

## [​](#serialization)Serialization

There are two key aspects to serialization in LangGraph:

{' ' * (self.list_depth - 1)}- `entrypoint` inputs and outputs must be JSON-serializable.

{' ' * (self.list_depth - 1)}- `task` outputs must be JSON-serializable.

These requirements are necessary for enabling checkpointing and workflow resumption. Use python primitives like dictionaries, lists, strings, numbers, and booleans to ensure that your inputs and outputs are serializable.
Serialization ensures that workflow state, such as task results and intermediate values, can be reliably saved and restored. This is critical for enabling human-in-the-loop interactions, fault tolerance, and parallel execution.
Providing non-serializable inputs or outputs will result in a runtime error when a workflow is configured with a checkpointer.

## [​](#determinism)Determinism

To utilize features like **human-in-the-loop**, any randomness should be encapsulated inside of **tasks**. This guarantees that when execution is halted (e.g., for human in the loop) and then resumed, it will follow the same *sequence of steps*, even if **task** results are non-deterministic.
LangGraph achieves this behavior by persisting **task** and [**subgraph**](/oss/python/langgraph/use-subgraphs) results as they execute. A well-designed workflow ensures that resuming execution follows the *same sequence of steps*, allowing previously computed results to be retrieved correctly without having to re-execute them. This is particularly useful for long-running **tasks** or **tasks** with non-deterministic results, as it avoids repeating previously done work and allows resuming from essentially the same.
While different runs of a workflow can produce different results, resuming a **specific** run should always follow the same sequence of recorded steps. This allows LangGraph to efficiently look up **task** and **subgraph** results that were executed prior to the graph being interrupted and avoid recomputing them.

## [​](#idempotency)Idempotency

Idempotency ensures that running the same operation multiple times produces the same result. This helps prevent duplicate API calls and redundant processing if a step is rerun due to a failure. Always place API calls inside **tasks** functions for checkpointing, and design them to be idempotent in case of re-execution. Re-execution can occur if a **task** starts, but does not complete successfully. Then, if the workflow is resumed, the **task** will run again. Use idempotency keys or verify existing results to avoid duplication.

## [​](#common-pitfalls)Common pitfalls

### [​](#handling-side-effects)Handling side effects

Encapsulate side effects (e.g., writing to a file, sending an email) in tasks to ensure they are not executed multiple times when resuming a workflow.

{' ' * (self.list_depth - 1)}- Incorrect
{' ' * (self.list_depth - 1)}- CorrectIn this example, a side effect (writing to a file) is directly included in the workflow, so it will be executed a second time when resuming the workflow.
```
@entrypoint(checkpointer=checkpointer)
def my_workflow(inputs: dict) -> int:
 # This code will be executed a second time when resuming the workflow.
 # Which is likely not what you want.
 with open("output.txt", "w") as f:
 f.write("Side effect executed")
 value = interrupt("question")
 return value

```
In this example, the side effect is encapsulated in a task, ensuring consistent execution upon resumption.
```
from langgraph.func import task

@task
def write_to_file():
 with open("output.txt", "w") as f:
 f.write("Side effect executed")

@entrypoint(checkpointer=checkpointer)
def my_workflow(inputs: dict) -> int:
 # The side effect is now encapsulated in a task.
 write_to_file().result()
 value = interrupt("question")
 return value

```

### [​](#non-deterministic-control-flow)Non-deterministic control flow

Operations that might give different results each time (like getting current time or random numbers) should be encapsulated in tasks to ensure that on resume, the same result is returned.

{' ' * (self.list_depth - 1)}- In a task: Get random number (5) → interrupt → resume → (returns 5 again) → …

{' ' * (self.list_depth - 1)}- Not in a task: Get random number (5) → interrupt → resume → get new random number (7) → …

This is especially important when using **human-in-the-loop** workflows with multiple interrupt calls. LangGraph keeps a list of resume values for each task/entrypoint. When an interrupt is encountered, it’s matched with the corresponding resume value. This matching is strictly **index-based**, so the order of the resume values should match the order of the interrupts.
If order of execution is not maintained when resuming, one [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) call may be matched with the wrong `resume` value, leading to incorrect results.
Please read the section on [determinism](#determinism) for more details.

{' ' * (self.list_depth - 1)}- Incorrect
{' ' * (self.list_depth - 1)}- CorrectIn this example, the workflow uses the current time to determine which task to execute. This is non-deterministic because the result of the workflow depends on the time at which it is executed.
```
from langgraph.func import entrypoint

@entrypoint(checkpointer=checkpointer)
def my_workflow(inputs: dict) -> int:
 t0 = inputs["t0"]
 t1 = time.time()

 delta_t = t1 - t0

 if delta_t > 1:
 result = slow_task(1).result()
 value = interrupt("question")
 else:
 result = slow_task(2).result()
 value = interrupt("question")

 return {
 "result": result,
 "value": value
 }

```
In this example, the workflow uses the input `t0` to determine which task to execute. This is deterministic because the result of the workflow depends only on the input.
```
import time

from langgraph.func import task

@task
def get_time() -> float:
 return time.time()

@entrypoint(checkpointer=checkpointer)
def my_workflow(inputs: dict) -> int:
 t0 = inputs["t0"]
 t1 = get_time().result()

 delta_t = t1 - t0

 if delta_t > 1:
 result = slow_task(1).result()
 value = interrupt("question")
 else:
 result = slow_task(2).result()
 value = interrupt("question")

 return {
 "result": result,
 "value": value
 }

```

## [​](#learn-more)Learn more

{' ' * (self.list_depth - 1)}- [How to use the Functional API](/oss/python/langgraph/use-functional-api)

{' ' * (self.list_depth - 1)}- [Graph API conceptual overview](/oss/python/langgraph/graph-api)

{' ' * (self.list_depth - 1)}- [Choosing between Graph API and Functional API](/oss/python/langgraph/choosing-apis)

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langgraph/functional-api.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Use the graph APIPrevious](/oss/python/langgraph/use-graph-api)[Use the functional APINext](/oss/python/langgraph/use-functional-api)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
